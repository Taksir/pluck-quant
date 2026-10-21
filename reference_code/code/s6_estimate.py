"""Stage 6: gene estimates per set and sample by ridge-regularised NNLS on the probe x transcript matrix.
Model: E[c_p] = sum_t n_pt * lam_t (n_pt = read starts in t whose 50 nt read covers probe p; lam_t = reads per start).
Gene estimate N_g = sum_{t in g} lam_t * max(L_t - R + 1, K). Tiny ridge (1e-6 x largest column norm^2) picks the
minimum-norm split where isoforms are not separable (for constitutive-only genes this is the current estimator)."""
import sys, pickle, json, time, numpy as np, scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.optimize import nnls
from numba import njit
@njit(cache=True)
def cdnnls_sp(indptr, indices, data, diag, Bm, L0, cols):
    # same coordinate descent, sparse symmetric Q in CSC form (large components)
    n = len(diag)
    for s in cols:
        x = np.maximum(L0[:, s], 0.0)
        g = Bm[:, s].copy()
        for i in range(n):
            if x[i] != 0.0:
                for k in range(indptr[i], indptr[i + 1]): g[indices[k]] -= data[k] * x[i]
        for it in range(100000):
            mc = 0.0; mx = 1e-12
            for i in range(n):
                new = max(0.0, x[i] + g[i] / diag[i])
                d = new - x[i]
                if d != 0.0:
                    for k in range(indptr[i], indptr[i + 1]): g[indices[k]] -= d * data[k]
                    x[i] = new
                    if abs(d) > mc: mc = abs(d)
                if x[i] > mx: mx = x[i]
            if mc <= 1e-10 * mx: break
        L0[:, s] = x
    return L0
@njit(cache=True)
def cdnnls(Q, Bm, L0, cols):
    # exact minimiser of the strictly convex (ridge) quadratic x'Qx - 2b'x over x >= 0, cyclic coordinate descent
    n = Q.shape[0]
    for s in cols:
        x = np.maximum(L0[:, s], 0.0)
        g = Bm[:, s] - Q @ x
        for it in range(20000):
            mc = 0.0; mx = 1e-12
            for i in range(n):
                new = max(0.0, x[i] + g[i] / Q[i, i])
                d = new - x[i]
                if d != 0.0:
                    g -= d * Q[:, i]; x[i] = new
                    if abs(d) > mc: mc = abs(d)
                if x[i] > mx: mx = x[i]
            if mc <= 1e-10 * mx: break
        L0[:, s] = x
    return L0
R, K = 50, 25
t0 = time.time()
D = np.load('w/design.npz')
import os; srrs = open(os.environ.get('SRRS', 'srrs.txt')).read().split()
Cn = np.stack([np.load(f'counts/{s}.npy').astype(np.float32) for s in srrs], 1)   # probes x samples
ltil = np.maximum(D['tx_len'] - R + 1, K).astype(np.float64)
ng = len(D['gene_id'])
sets = sorted({k[4:].rsplit('__', 1)[0] for k in D.files if k.startswith('set_')})
if os.environ.get('ONLY'): sets = os.environ['ONLY'].split(',')
out = {}
errs = []
for name in sets:
    T1 = name.startswith('T1') or name == 'base'
    txok = D['tx_tsl1pc'] if T1 else np.ones(len(D['tx_len']), bool)
    for v in ('ref', 'vf'):
        pr = D[f'set_{name}__{v}']
        inset = np.zeros(len(D['sense_code']), bool); inset[pr] = True
        m = inset[D['A_probe']] & txok[D['A_tx']] & (D['A_nstarts'] > 0)
        rows, cols, val = D['A_probe'][m], D['A_tx'][m], D['A_nstarts'][m].astype(float)
        up, ri = np.unique(rows, return_inverse=True); ut, ci = np.unique(cols, return_inverse=True)
        M = sp.csr_matrix((val, (ri, ci)), shape=(len(up), len(ut)))
        nc, lab = connected_components((M.T @ M), directed=False)
        Mc = M.tocsc()
        lam = np.zeros((len(ut), len(srrs)))
        order = np.argsort(lab, kind='stable'); bnd = np.r_[0, np.cumsum(np.bincount(lab, minlength=nc))]
        Y = Cn[up]
        n_nnls = 0; maxcomp = 0
        for c in range(nc):
            tc = order[bnd[c]:bnd[c + 1]]
            sub = Mc[:, tc]
            rr = np.unique(sub.indices)
            y = Y[rr].astype(np.float64)
            maxcomp = max(maxcomp, len(tc))
            if len(tc) > 300:
                Xs = sub[rr]
                XtXs = (Xs.T @ Xs).tocsc(); a = 1e-6 * XtXs.diagonal().max()
                Qs = (XtXs + a * sp.identity(len(tc), format='csc')).tocsc(); Bm = np.asarray(Xs.T @ y)
                from scipy.sparse.linalg import splu
                L = splu(Qs).solve(Bm)
                bad = np.nonzero((L < 0).any(0))[0]
                if len(bad):
                    Qs.sort_indices()
                    L = cdnnls_sp(Qs.indptr, Qs.indices, Qs.data, Qs.diagonal(), Bm, np.ascontiguousarray(L), bad)
                    n_nnls += len(bad)
                lam[tc] = L
                continue
            X = sub[rr].toarray()
            XtX = X.T @ X; a = 1e-6 * XtX.diagonal().max()
            L = np.linalg.solve(XtX + a * np.eye(len(tc)), X.T @ y)
            bad = np.nonzero((L < 0).any(0))[0]
            if len(bad):
                Q = XtX + a * np.eye(len(tc)); Bm = X.T @ y
                L = cdnnls(Q, Bm, np.ascontiguousarray(L), bad)
                if len(errs) < 300 and c % 50 == 0:   # spot-check against scipy's Lawson-Hanson
                    Xa = np.vstack([X, np.sqrt(a) * np.eye(len(tc))])
                    ref = nnls(Xa, np.r_[y[:, bad[0]], np.zeros(len(tc))], maxiter=50 * len(tc))[0]
                    gref = (ref * ltil[ut[tc]]).sum(); gcd = (L[:, bad[0]] * ltil[ut[tc]]).sum()
                    errs.append(abs(gref - gcd) / max(gref, 1e-9))
                n_nnls += len(bad)
            lam[tc] = L
        gene_of_t = D['tx_gene'][ut]
        est = np.zeros((ng, len(srrs)))
        np.add.at(est, gene_of_t, lam * ltil[ut][:, None])
        has = np.zeros(ng, bool); has[gene_of_t] = True
        # gene-specific probe: all its compatible model transcripts belong to one gene
        pg = gene_of_t[ci]
        first = np.full(len(up), -1); np.maximum.at(first, ri, pg); mn = np.full(len(up), 1 << 30); np.minimum.at(mn, ri, pg)
        spec_g = np.zeros(ng, bool); spec_g[first[first == mn]] = True
        out[(name, v)] = dict(est=est, has=has, specific=spec_g, n_probes=len(up), n_tx=len(ut), n_comp=nc, max_comp=maxcomp, n_nnls=n_nnls)
        print(name, v, 'probes', len(up), 'tx', len(ut), 'genes', int(has.sum()), 'comps', nc, 'max', maxcomp, 'nnls', n_nnls, round(time.time() - t0), flush=True)
pickle.dump(dict(srrs=srrs, out=out), open(os.environ.get('EST', 'w/estimates.pkl'), 'wb'))
print('S6DONE component-total rel diff CD vs scipy nnls: max', max(errs) if errs else None, 'n', len(errs), flush=True)
