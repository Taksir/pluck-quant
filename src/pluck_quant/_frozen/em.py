"""Sparse-probe EM over transcripts with an optional per-sample positional (5'/3') and GC bias model.
Counts c_k ~ Poisson(sum_t a_kt X_t); a_kt = w_s(pos, L_t, gc_k) * (R-K+1) / Z_t, Z_t = sum over all k-mer starts of w."""
import numpy as np, numba
K, R = 25, 50
NB = 24  # distance bins (25 nt) for 5' and 3' profiles; last bin = >= 575 nt
def dbin(d):
    return np.minimum(d // 25, NB - 1).astype(np.int64)

@numba.njit(cache=True)
def em_run(c, kp, kt, ka, ntx, A, niter, tol):
    """kp/kt/ka: compat entries (probe, transcript, rate), sorted by probe. A[t] = sum of rates of t over used probes."""
    X = np.ones(ntx)
    n = len(kp)
    for it in range(niter):
        acc = np.zeros(ntx)
        i = 0
        while i < n:
            k = kp[i]; j = i; mu = 0.0
            while j < n and kp[j] == k:
                mu += ka[j] * X[kt[j]]; j += 1
            if mu > 0 and c[k] > 0:
                f = c[k] / mu
                for q in range(i, j):
                    acc[kt[q]] += ka[q] * X[kt[q]] * f
            i = j
        delta = 0.0
        for t in range(ntx):
            if A[t] > 0:
                nx = acc[t] / A[t]
                if nx > 1.0:
                    d = abs(nx - X[t]) / nx
                    if d > delta: delta = d
                X[t] = nx
            else:
                X[t] = 0.0
        if delta < tol:
            break
    return X, it

class Model:
    def __init__(self, D, use, gc):
        """D: design npz dict; use: bool over union probes allowed in EM; gc: GC count (0..25) per union probe."""
        m = use[D['cp']]
        self.kp = D['cp'][m]; self.kt = D['ct'][m]; self.pos = D['cpos'][m].astype(np.int64)
        o = np.argsort(self.kp, kind='stable')
        self.kp, self.kt, self.pos = self.kp[o], self.kt[o], self.pos[o]
        self.L = D['tlen'].astype(np.int64); self.ntx = len(self.L)
        self.d5 = dbin(self.pos); self.d3 = dbin(self.L[self.kt] - self.pos - K)
        self.gcb = np.minimum(gc[self.kp] // 3, 7)
        # single-transcript probes (for bias profiles)
        nt = np.bincount(self.kp, minlength=len(use))
        self.single = nt[self.kp] == 1
        # per transcript, number of k-mer starts in each (d5, d3) bin pair -> for Z_t
        self.Lk = np.maximum(self.L - K + 1, 0)
        self.prof = None
        # equivalence classes: probes with identical (transcript, d5, d3, gc-bin) entry lists
        code = ((self.kt.astype(np.int64) * NB + self.d5) * NB + self.d3) * 8 + self.gcb
        start = np.searchsorted(self.kp, np.arange(len(use) + 1))
        h1, h2 = phash(code, start)
        has = start[1:] > start[:-1]
        key = np.stack([h1, h2], 1)
        _, first, ec = np.unique(key, axis=0, return_index=True, return_inverse=True)
        ec = ec.ravel()
        self.ec_of_probe = np.where(has, ec, -1)
        nec = len(first)
        rep = first  # representative probe per EC
        cnt = start[rep + 1] - start[rep]
        ent = np.concatenate([np.arange(start[p], start[p + 1]) for p in rep]) if nec < 0 else rep_entries(start, rep)
        self.e_ec = np.repeat(np.arange(nec), cnt); self.e_ix = ent  # EC entry -> probe-entry index (for rates)
        self.nec = nec
    def rates(self, h5, h3, hg):
        w = h5[self.d5] * h3[self.d3] * hg[self.gcb]
        Z = zsum(self.Lk, h5, h3)
        Z = np.maximum(Z, 1.0)
        a = w * (R - K + 1) / Z[self.kt]
        return a
    def fit(self, c, bias=True, niter=3000, tol=1e-3, rounds=2, frag=None, X0=None):
        h5 = np.ones(NB); h3 = np.ones(NB); hg = np.ones(8)
        if frag is not None:  # analytic 5' edge from a normal fragment-length distribution (mean, sd), Salmon's SE default 250/25
            from scipy.stats import norm
            mu, sd = frag
            ctr = np.arange(NB) * 25 + 12
            h5 = np.array([norm.cdf((p + np.arange(K, R + 1) - mu) / sd).mean() for p in ctr]); h5 = np.maximum(h5, 1e-3)
        m = self.ec_of_probe >= 0
        cec = np.bincount(self.ec_of_probe[m], weights=c[m].astype(np.float64), minlength=self.nec)
        nprb = np.bincount(self.ec_of_probe[m], minlength=self.nec).astype(np.float64)
        for r in range(rounds if bias else 1):
            a = self.rates(h5, h3, hg)
            A = np.bincount(self.kt, weights=a, minlength=self.ntx)
            ae = a[self.e_ix] * nprb[self.e_ec]
            X, it = squarem(cec, self.e_ec, self.kt[self.e_ix], ae, self.ntx, A, niter, tol, X0=(X if r > 0 else X0))
            if not bias or r == rounds - 1:
                break
            # update profiles from single-transcript probes: obs / expected (with current w removed), alternating
            s = self.single & (X[self.kt] > 20)
            base = (R - K + 1) / np.maximum(zsum(self.Lk, h5, h3), 1.0)[self.kt[s]] * X[self.kt[s]]
            obs = c[self.kp[s]].astype(np.float64)
            d5, d3, gb = self.d5[s], self.d3[s], self.gcb[s]
            for _ in range(5):
                e = base * h3[d3] * hg[gb]; h5 = np.bincount(d5, obs, NB) / np.maximum(np.bincount(d5, e, NB), 1e-9)
                h5 /= h5[-1] if h5[-1] > 0 else 1
                e = base * h5[d5] * hg[gb]; h3 = np.bincount(d3, obs, NB) / np.maximum(np.bincount(d3, e, NB), 1e-9)
                h3 /= h3[-1] if h3[-1] > 0 else 1
                e = base * h5[d5] * h3[d3]; hg = np.bincount(gb, obs, 8) / np.maximum(np.bincount(gb, e, 8), 1e-9)
                hg /= np.average(hg, weights=np.bincount(gb, minlength=8) + 1e-9)
            h5 = np.clip(h5, 0.02, 50); h3 = np.clip(h3, 0.02, 50); hg = np.clip(hg, 0.05, 20)
        self.prof = (h5, h3, hg)
        return X, it

@numba.njit(cache=True)
def zsum(Lk, h5, h3):
    nb = len(h5)
    Z = np.zeros(len(Lk))
    for t in range(len(Lk)):
        L = Lk[t]; s = 0.0
        for p in range(L):
            b5 = p // 25
            if b5 > nb - 1: b5 = nb - 1
            b3 = (L - 1 - p) // 25
            if b3 > nb - 1: b3 = nb - 1
            s += h5[b5] * h3[b3]
        Z[t] = s
    return Z

@numba.njit(cache=True)
def phash(code, start):
    n = len(start) - 1
    h1 = np.zeros(n, np.int64); h2 = np.zeros(n, np.int64)
    for p in range(n):
        a = np.int64(1469598103934665603); b = np.int64(7)
        for i in range(start[p], start[p + 1]):
            a = (a ^ code[i]) * np.int64(1099511628211)
            b = b * np.int64(1000003) + code[i] * np.int64(31) + np.int64(17)
        h1[p] = a; h2[p] = b
    return h1, h2

@numba.njit(cache=True)
def rep_entries(start, rep):
    tot = 0
    for p in rep: tot += start[p + 1] - start[p]
    out = np.empty(tot, np.int64); j = 0
    for p in rep:
        for i in range(start[p], start[p + 1]):
            out[j] = i; j += 1
    return out

@numba.njit(cache=True)
def em_step(c, kp, kt, ka, ntx, A, X):
    acc = np.zeros(ntx)
    n = len(kp); i = 0
    while i < n:
        k = kp[i]; j = i; mu = 0.0
        while j < n and kp[j] == k:
            mu += ka[j] * X[kt[j]]; j += 1
        if mu > 0 and c[k] > 0:
            f = c[k] / mu
            for q in range(i, j):
                acc[kt[q]] += ka[q] * X[kt[q]] * f
        i = j
    out = np.zeros(ntx)
    for t in range(ntx):
        if A[t] > 0:
            out[t] = acc[t] / A[t]
    return out

GENE_OF_TX = None  # set by caller: int gene index per transcript; convergence is checked on gene sums
def squarem(c, kp, kt, ka, ntx, A, niter, tol, X0=None):
    """SQUAREM-accelerated EM (Varadhan & Roland 2008); stop when the max relative change of gene sums > 10 is < tol."""
    X = np.where(A > 0, 1.0, 0.0) if X0 is None else np.where(A > 0, np.maximum(X0, 1e-3), 0.0)
    gsum = lambda x: np.bincount(GENE_OF_TX, x) if GENE_OF_TX is not None else x
    prev = gsum(X)
    for it in range(niter):
        X1 = em_step(c, kp, kt, ka, ntx, A, X)
        X2 = em_step(c, kp, kt, ka, ntx, A, X1)
        r = X1 - X; v = X2 - X1 - r
        nv = np.sqrt((v * v).sum())
        a = -np.sqrt((r * r).sum()) / nv if nv > 0 else -1.0
        a = min(a, -1.0)
        Xn = np.maximum(X - 2 * a * r + a * a * v, 0.0)
        Xn = em_step(c, kp, kt, ka, ntx, A, Xn)
        g = gsum(Xn); big = g > 10
        d = np.max(np.abs(g[big] - prev[big]) / g[big]) if big.any() else 0.0
        X = Xn; prev = g
        if d < tol:
            break
    return X, it
