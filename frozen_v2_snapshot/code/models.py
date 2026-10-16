# Gene-probe model families (26 nt vf probes) + EM variants, strain-grouped 5-fold held-out evaluation on Keele liver.
import json, os, sys, numpy as np, pandas as pd
from common import *
import score as SC
D = np.load('design.npz'); gc = np.load('gc.npy')
m = meta(); srrs = sorted(m.srr.unique()); ns = len(srrs)
rng = np.random.default_rng(20261005); fold = np.zeros(ns, int); fold[rng.permutation(ns)] = np.arange(ns) % 5
C = np.stack([np.load(f'counts/{s}.npy') for s in srrs]).astype(np.float64)
vf = D['A_vf']; ai = D['A_idx'][vf]; ag = D['A_gene'][vf].astype(np.int64); apos = D['A_pos'][vf].astype(np.int64)
pc = C[:, ai]; del C
ng = len(G); mprobe = np.bincount(ag, minlength=ng)
scale = np.where(mprobe > 0, LT / (np.maximum(mprobe, 1) * (R - K + 1)), np.nan)
S = np.stack([np.bincount(ag, weights=pc[i], minlength=ng) for i in range(ns)])
relb = np.minimum((apos / LREF[ag] * 10).astype(int), 9)
d5b = np.minimum(apos // 50, 9); d3b = np.minimum(np.maximum(LREF[ag] - apos - K, 0) // 50, 9)
gcb = np.minimum(gc[ai] // 3, 7)
EST = {'vf26': S * scale}

def rho_from(train):
    tot = pc[train].sum(0); St = S[train].sum(0)
    return np.where(St[ag] > 0, tot / np.maximum(St[ag], 1e-9), 1.0 / np.maximum(mprobe[ag], 1))

def robust(i, rho, f=None, z=5.0, iters=3):
    c = pc[i]; w = rho * (f if f is not None else 1.0)
    keep = np.ones(len(c), bool)
    th = np.bincount(ag, c, ng) / np.maximum(np.bincount(ag, w, ng), 1e-12)
    for _ in range(iters):
        e = w * th[ag]
        bad = ((c == 0) & (e >= 10)) | (np.abs(c - e) > z * np.sqrt(np.maximum(e, 1)) * 2)  # 2x allows modest overdispersion
        keep = ~bad
        th = np.bincount(ag, c * keep, ng) / np.maximum(np.bincount(ag, w * keep, ng), 1e-12)
    return th, 1 - keep.mean()

def shape(i, rho, th):
    """per-sample multiplicative factors over relative position decile x GC bin, from genes with theta >= 50"""
    c = pc[i]; e0 = rho * th[ag]; ok = th[ag] >= 50
    fp = np.ones(10); fg = np.ones(8)
    for _ in range(6):
        e = e0 * fg[gcb]; fp = np.bincount(relb[ok], c[ok], 10) / np.maximum(np.bincount(relb[ok], e[ok], 10), 1e-9)
        e = e0 * fp[relb]; fg = np.bincount(gcb[ok], c[ok], 8) / np.maximum(np.bincount(gcb[ok], e[ok], 8), 1e-9)
    return fp[relb] * fg[gcb], fp, fg

out_est = {k: np.full((ns, ng), np.nan) for k in ('robust', 'shape', 'shape_robust', 'shrink')}
drop = []; shapes = {}
for k in range(5):
    tr = np.nonzero(fold != k)[0]; te = np.nonzero(fold == k)[0]
    rho = rho_from(tr)
    # shrinkage parameters from training strains: per gene mean and between-strain variance of depth-normalised log estimate
    N = S.sum(1)
    y = np.log(np.maximum(S, 0.5)) - np.log(N)[:, None]
    mu = y[tr].mean(0); s2 = (1 / np.maximum(S[tr], 1)).mean(0); tau2 = np.maximum(y[tr].var(0) - s2, 1e-4)
    for i in te:
        th, fr = robust(i, rho); out_est['robust'][i] = th * scale; drop.append(fr)
        th0 = S[i]
        f, fp, fg = shape(i, rho, th0); shapes[srrs[i]] = dict(pos=fp.round(3).tolist(), gc=fg.round(3).tolist())
        out_est['shape'][i] = np.bincount(ag, pc[i], ng) / np.maximum(np.bincount(ag, rho * f, ng), 1e-12) * scale
        th2, _ = robust(i, rho, f); out_est['shape_robust'][i] = th2 * scale
        B = tau2 / (tau2 + 1 / np.maximum(S[i], 1))
        ys = mu + B * (y[i] - mu)
        out_est['shrink'][i] = np.exp(ys + np.log(N[i])) * scale
EST.update(out_est)
# EM variants
for v in ('emu', 'emf', 'eme', 'emt', 'emg'):
    if all(os.path.exists(f'emest/{s}_{v}.npz') for s in srrs):
        EST[v] = np.stack([np.load(f'emest/{s}_{v}.npz')['gene'] for s in srrs])
# Hybrid: gene probes where their nominal read capture is high, EM elsewhere. Threshold tau picked on training folds
cap = np.where(mprobe > 0, mprobe * (R - K + 1) / LT, 0)
mad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v, 0)), 0)
def wg_sd(est, rows):
    # within-gene robust SD of calibrated log error on training rows (truth >= 50), median over genes with >= 10 rows
    tr = np.full((ns, ng), np.nan); tv = m.pivot_table(index='srr', columns='gi', values='truth').reindex(srrs)
    tr[:, tv.columns.to_numpy()] = tv.to_numpy()
    e = est[rows]; t = tr[rows]
    Cs = np.nanmedian(np.where((t >= 50) & (e > 0), t / np.where(e > 0, e, 1), np.nan), 1)
    lr = np.where(t >= 50, np.log2((Cs[:, None] * e + 1) / (t + 1)), np.nan)
    ok = (np.isfinite(lr).sum(0) >= 10)
    return float(np.nanmedian(mad(lr[:, ok])))
TAUS = [0.0, 0.3, 0.45, 0.6, 0.75, 0.9, 99]
hyb_info = {}
for em_name in [v for v in ('emu', 'emt', 'emg') if v in EST]:
    h = np.full((ns, ng), np.nan); picked = []
    for k in range(5):
        tr = np.nonzero(fold != k)[0]; te = np.nonzero(fold == k)[0]
        sc = {t: wg_sd(np.where(cap[None, :] >= t, EST['vf26'], EST[em_name]), tr) for t in TAUS}
        tb = min(sc, key=sc.get); picked.append(tb)
        h[te] = np.where(cap[None, :] >= tb, EST['vf26'][te], EST[em_name][te])
    EST['hyb_' + em_name] = h; hyb_info[em_name] = picked
    # observed capture: share of the gene's EM-estimated fragments that its gene-unique constitutive probes see (training strains, no truth)
    h2 = np.full((ns, ng), np.nan); picked2 = []
    for k in range(5):
        tr = np.nonzero(fold != k)[0]; te = np.nonzero(fold == k)[0]
        fo = np.nansum(S[tr], 0) / np.maximum(np.nansum(EST[em_name][tr], 0), 1e-9)
        sc = {t: wg_sd(np.where(fo[None, :] >= t, EST['vf26'], EST[em_name]), tr) for t in TAUS}
        tb = min(sc, key=sc.get); picked2.append(tb)
        h2[te] = np.where(fo[None, :] >= tb, EST['vf26'][te], EST[em_name][te])
    EST['hyb2_' + em_name] = h2; hyb_info[em_name + '_observed'] = picked2
# Salmon-trained per-gene calibration (absolute level only), fold-held-out, for vf26 and best EM
sal = m.pivot_table(index='srr', columns='gi', values='salmon_pan').reindex(srrs)
salm = np.full((ns, ng), np.nan); salm[:, sal.columns.to_numpy()] = sal.to_numpy()
for base in [b for b in ('vf26', 'emu', 'emg', 'hyb_emu', 'hyb_emg') if b in EST]:
    cal = np.full((ns, ng), np.nan)
    for k in range(5):
        tr = fold != k; te = fold == k
        with np.errstate(all='ignore'):
            lr = np.log((salm[tr] + 1) / (EST[base][tr] + 1))
        kg = np.nanmedian(lr, 0); kg = np.where(np.isfinite(kg), kg, 0)
        cal[te] = (EST[base][te] + 1) * np.exp(kg) - 1
    EST[base + '_salcal'] = cal
for k, v in EST.items():
    m = gene_matrix(m, v, k)
meths = list(EST) + ['salmon_ref', 'salmon_pan']
res, d = SC.score(m, meths)
res['_info'] = dict(robust_drop_frac=float(np.mean(drop)), folds=fold.tolist(), hybrid_tau_per_fold=hyb_info)
np.savez_compressed('estimates.npz', srrs=np.array(srrs), **EST)
json.dump(res, open('models.json', 'w'), indent=1)
json.dump(shapes, open('shapes.json', 'w'))
d.to_parquet('scored.parquet') if False else None
pd.DataFrame(res).T.to_csv('models.tsv', sep='\t')
print(pd.DataFrame({k: v for k, v in res.items() if not k.startswith('_')}).T.to_string())
