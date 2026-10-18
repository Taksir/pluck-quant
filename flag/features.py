"""Design-only per-gene features for the frozen core tier (no sample, no truth, no Salmon).
usage: features.py 26|76 OUT.tsv.gz
  26 -> probe_em/design.npz (+ out/genome_count.npy), gene probes probes_s26, R = 50  (Keele, 50 bp SE)
  76 -> e5_munger/design/dsn76 (EM set at SP = 76), gene probes probes_s76 (R = 100, Munger) and probes_s77 (R = 101, GSE199702)
EM set = frozen emg rule (run_emg.py): inB & genome_copies <= max(1, n_genes_containing) & n_genes_containing == 1.
Gene-probe incidence in non-TSL1 transcripts is taken from the 26 nt design (gene probes are constitutive in TSL1 isoforms
at every spacing, so the per-transcript fraction m_t/m is a property of exon structure)."""
import sys, json, numpy as np, pandas as pd
from scipy.optimize import linprog
R_ = '/mnt/project-files/revision'; K = 25
MODE, OUT = sys.argv[1], sys.argv[2]
G = json.load(open(f'{R_}/fast_competitors/e9_prep_cache/ref/genes.json'))
GID = np.array([g['gene_id'] for g in G]); ng = len(G)
gpos = {g: i for i, g in enumerate(GID)}
if MODE == '26':
    D = dict(np.load(f'{R_}/probe_em/design.npz')); gcount = np.load(f'{R_}/probe_em/out/genome_count.npy')
    EMS, EMR = 26, 50
    GP = {'26': (f'{R_}/respacing/probes/probes_s26.npz', 50)}
else:
    D = dict(np.load(f'{R_}/e5_munger/design/dsn76/design.npz')); gcount = np.load(f'{R_}/e5_munger/design/dsn76/genome_count.npy')
    EMS, EMR = 76, 100
    GP = {'76': (f'{R_}/respacing/../e5_munger/probes_s76.npz', 100), '77': (f'{R_}/respacing/probes/probes_s77.npz', 101)}
tg_u, tg_i = np.unique(D['tgene'], return_inverse=True); tg_i = tg_i.ravel()
gmap = np.array([gpos.get(g, -1) for g in tg_u])          # tgene index -> genes.json index
txg = gmap[tg_i]                                          # transcript -> genes.json index (-1 outside core universe)
cp, ct = D['cp'], D['ct']; nU = len(D['U'])
pairs = np.unique(np.stack([cp, tg_i[ct]], 1), axis=0); ngp = np.bincount(pairs[:, 0], minlength=nU)
inB = D['inB']; use = inB & (gcount <= np.maximum(ngp, 1)); emg = use & (ngp == 1)
L = D['tlen'].astype(np.int64); bio = D['tbio']; ntxall = len(L)
F = pd.DataFrame({'gene_id': GID})
# ---- transcript structure
gi_tx = txg >= 0
F['ntx'] = np.bincount(txg[gi_tx], minlength=ng)
F['ntx_pc'] = np.bincount(txg[gi_tx & (bio == 'protein_coding')], minlength=ng)
F['ntx_ri_nmd'] = np.bincount(txg[gi_tx & np.isin(bio, ['retained_intron', 'nonsense_mediated_decay', 'processed_transcript'])], minlength=ng)
F['ntsl1'] = [len(g['tx_lens']) for g in G]
F['tsl1_len_logratio'] = [np.log2(max(g['tx_lens']) / min(g['tx_lens'])) for g in G]
F['ref_len'] = [g['tx_lens'][0] for g in G]
# ---- B windows (all variant-free candidate probes of the gene's transcripts): shared / multi-copy content
eB = inB[cp] & gi_tx[ct]
pg = np.unique(np.stack([cp[eB], txg[ct[eB]]], 1), axis=0)   # (probe, gene) unique
p, g = pg[:, 0], pg[:, 1]
F['nB'] = np.bincount(g, minlength=ng)
F['nB_shared'] = np.bincount(g, weights=(ngp[p] > 1), minlength=ng)
F['nB_multicopy'] = np.bincount(g, weights=(gcount[p] > np.maximum(ngp[p], 1)), minlength=ng)
F['n_emg'] = np.bincount(g, weights=emg[p], minlength=ng)
with np.errstate(invalid='ignore', divide='ignore'):
    F['frac_shared'] = F.nB_shared / F.nB; F['frac_multicopy'] = F.nB_multicopy / F.nB; F['frac_emg'] = F.n_emg / F.nB
# ---- per-transcript EM probe occurrences and EM expected capture (uniform reads per transcript)
eE = emg[cp] & gi_tx[ct]
nt = np.bincount(ct[eE], minlength=ntxall)              # emg probe occurrences per transcript
Leff_em = np.maximum(L - EMR + 1, K).astype(float)
cap_t = np.minimum(nt * EMS / Leff_em, 1.0)
F['tx_invisible'] = np.bincount(txg[gi_tx], weights=(nt[gi_tx] == 0), minlength=ng)
F['frac_tx_invisible'] = F.tx_invisible / F.ntx.clip(lower=1)
pcm = gi_tx & (bio == 'protein_coding')
F['frac_pc_invisible'] = np.bincount(txg[pcm], weights=(nt[pcm] == 0), minlength=ng) / F.ntx_pc.clip(lower=1)
F['em_capture'] = np.bincount(txg[gi_tx], weights=cap_t[gi_tx], minlength=ng) / F.ntx.clip(lower=1)
# ---- gene probes: counts, variant loss, expected capture, predicted observed capture f (hybrid switch)
D26 = np.load(f'{R_}/probe_em/design.npz') if MODE != '26' else None
A = D if MODE == '26' else {k: D26[k] for k in ('A_idx', 'A_gene', 'A_vf', 'cp', 'ct', 'U')}
# m_t / m from the 26 nt design: fraction of the gene's vf gene probes present in each transcript
Avf = A['A_vf']; aidx = A['A_idx'][Avf]; agene = A['A_gene'][Avf].astype(np.int64)
isA = np.zeros(len(A['U']), np.int64) - 1; isA[aidx] = agene
eA = isA[A['cp']] >= 0
tA = A['ct'][eA]; gA = isA[A['cp'][eA]]
own = txg[tA] == gA                                     # gene probe found in a transcript of its own gene
mt = np.bincount(tA[own], minlength=ntxall).astype(float)
m26 = np.bincount(agene, minlength=ng).astype(float)
for sp, (fn, RR) in GP.items():
    Z = np.load(fn); gg = Z['gene'].astype(np.int64)
    mvf = np.bincount(gg[Z['in_vf']], minlength=ng); mref = np.bincount(gg[Z['in_ref']], minlength=ng)
    S = RR - K + 1
    lt = np.array([np.mean([max(x - RR + 1, K) for x in g['tx_lens']]) for g in G])
    F[f'm_vf_{sp}'] = mvf; F[f'm_ref_{sp}'] = mref
    F[f'vf_loss_{sp}'] = np.where(mref > 0, 1 - mvf / np.maximum(mref, 1), np.nan)
    F[f'gp_capture_{sp}'] = np.minimum(mvf * S / lt, 1.0)
    Le = np.maximum(L - RR + 1, K).astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        contrib = (mt / np.maximum(m26[np.maximum(txg, 0)], 1)) * lt[np.maximum(txg, 0)] / Le
    vis = gi_tx & (nt > 0)                               # transcripts the EM can see
    for nm, msk in (('all', vis), ('pc', vis & (bio == 'protein_coding'))):
        num = np.bincount(txg[msk], weights=contrib[msk], minlength=ng); den = np.bincount(txg[msk], minlength=ng)
        F[f'f_pred_{nm}_{sp}'] = np.where((den > 0) & (m26 > 0), num / np.maximum(den, 1), np.nan)
# ---- EM identifiability per gene: probe patterns over the gene's transcripts
# rows = distinct transcript-occurrence patterns of emg probes (weight = number of probes), entries occ * S / Leff
o = np.argsort(cp[eE], kind='stable'); P = cp[eE][o]; T = ct[eE][o]
rng = np.random.default_rng(1); hv = rng.integers(1, 2**62, ntxall, dtype=np.int64)
st = np.r_[0, np.nonzero(np.diff(P))[0] + 1]
sig = np.add.reduceat(hv[T], st)                        # order-free pattern signature (multiplicity-aware)
pg2 = txg[T[st]]
key = pd.DataFrame({'g': pg2, 'sig': sig, 'st': st})
grp = key.groupby(['g', 'sig'], sort=True).agg(st=('st', 'first'), w=('st', 'size')).reset_index()
en = np.r_[st[1:], len(P)]; end_of = dict(zip(st.tolist(), en.tolist()))
amb = np.full(ng, np.nan); vfac = np.full(ng, np.nan); nullfr = np.full(ng, np.nan); cond = np.full(ng, np.nan); nvis = np.zeros(ng, int); rank_def = np.zeros(ng, int)
for gi, sub in grp.groupby('g', sort=False):
    txs = np.unique(np.concatenate([T[s:end_of[s]] for s in sub.st]))
    ix = {t: j for j, t in enumerate(txs)}; nv = len(txs); nvis[gi] = nv
    rows = np.zeros((len(sub), nv))
    for r, s in enumerate(sub.st):
        for t in T[s:end_of[s]]: rows[r, ix[t]] += EMS / Leff_em[t]
    w = sub.w.to_numpy(float)
    X0 = np.ones(nv); mu = rows @ X0
    I = (rows * (w / mu)[:, None]).T @ rows
    ev, V = np.linalg.eigh(I); tol = 1e-9 * max(ev.max(), 1e-300); good = ev > tol
    rk = int(good.sum()); rank_def[gi] = nv - rk
    pr = (V.T @ np.ones(nv)) ** 2
    vfac[gi] = float((pr[good] / ev[good]).sum()) / nv   # Var(gene sum) / Poisson(all reads), identifiable part
    nullfr[gi] = float(pr[~good].sum()) / nv              # share of the gene-sum direction in the unidentified space
    cond[gi] = float(np.log10(ev.max() / max(ev[good].min(), 1e-300)))
    if rk == nv:
        amb[gi] = 0.0
    else:
        b = mu
        lo = linprog(np.ones(nv), A_eq=rows, b_eq=b, bounds=(0, None), method='highs')
        hi = linprog(-np.ones(nv), A_eq=rows, b_eq=b, bounds=(0, None), method='highs')
        if lo.status == 0 and hi.status == 0: amb[gi] = (-hi.fun - lo.fun) / nv
        else: amb[gi] = np.inf if hi.status == 3 else np.nan
F['em_nvis'] = nvis; F['em_rank_def'] = rank_def; F['em_var_factor'] = vfac; F['em_ambiguity'] = amb
F['em_npat'] = np.bincount(grp.g, minlength=ng); F['em_null_frac'] = nullfr; F['em_log10_cond'] = cond
F.insert(1, 'design', MODE)
F.to_csv(OUT, sep='\t', index=False, float_format='%.5g')
print(OUT, F.shape, flush=True)
print(F.describe().T.to_string())
