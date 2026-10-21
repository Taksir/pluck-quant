#!/usr/bin/env python3
"""Gene-bootstrap 95% intervals for the simulated per-founder median shift (Fig. 2b), which the stored table gives as medians only.

Reads (read only):  revision/e2_sim/out/est/clean_{se50,se101,pe}.npz, revision/e2_sim/work/design.npz,
                    revision/respacing/probes/probes_s{26,77}.npz, revision/fast_competitors/e9_prep_cache/ref/genes.json
                    revision/e2_sim/out/delta_founder.tsv (for the reproduction check)
Writes:             figures/data/fig2b_e2_delta_founder_ci.tsv

Per-gene shift follows score_e2.py (calib, hit_matrix, m3) unchanged: for each hit gene of founder k (a reference-unique probe of the
layout's spacing overlaps a variant of k), delta = mean calibrated log2 error over samples carrying k minus mean over B6 samples
(at least 1 and 2 samples, truth >= 50). The median over genes must equal delta_founder.tsv to 4 dp; the script stops otherwise.
Interval: percentile bootstrap over genes, 2,000 resamples, seed 20261015.
"""
import json, sys
import numpy as np, pandas as pd

REV = '/mnt/project-files/revision'
OUT = '/mnt/project-files/paper/figures/data/fig2b_e2_delta_founder_ci.tsv'
NS, MIN_TRUTH, LET = 8, 50, 'ACDEFGH'
FBIT = {l: 1 << i for i, l in enumerate(LET)}
D = np.load(f'{REV}/e2_sim/work/design.npz')
SIMG = D['genes'].tolist(); sg = {g: i for i, g in enumerate(SIMG)}; NG = len(SIMG)
G = json.load(open(f'{REV}/fast_competitors/e9_prep_cache/ref/genes.json'))
J2S = np.array([sg.get(g['gene_id'], -1) for g in G])
hap = D['hap'].astype(int); B6 = hap == 0
T = np.stack([np.bincount(D['gene_idx'], weights=D['expected'][s], minlength=NG) for s in range(NS)])
stored = pd.read_csv(f'{REV}/e2_sim/out/delta_founder.tsv', sep='\t')
rng = np.random.default_rng(20261015)


def hit_matrix(R):
    P = np.load(f'{REV}/respacing/probes/probes_s{26 if R == 50 else 77}.npz')
    ng = len(G); h = np.zeros((ng, 7), bool); sel = P['in_ref']
    for i, l in enumerate(LET):
        h[:, i] = np.bincount(P['gene'][sel], weights=((P['mask'][sel] & FBIT[l]) > 0), minlength=ng) > 0
    out = np.zeros((NG, 7), bool); ok = J2S >= 0; out[J2S[ok]] = h[ok]
    return out


def calib(E, gset):
    C = np.full(NS, np.nan)
    for s in range(NS):
        q = B6[s] & gset & (T[s] >= MIN_TRUTH) & (E[s] > 0) & np.isfinite(E[s])
        C[s] = np.median(T[s, q] / E[s, q])
    return np.log2((C[:, None] * E + 1) / (T + 1))


rows = []
for lay, R in (('se50', 50), ('se101', 101), ('pe', 101)):
    E = dict(np.load(f'{REV}/e2_sim/out/est/clean_{lay}.npz'))
    HIT = hit_matrix(R)
    core = np.all(np.isfinite(E['hybrid']), 0)
    for m in ('hybrid', 'vf', 'emg', 'ref', 'salmon_ref', 'salmon_pan'):
        gset = core & np.all(np.isfinite(E[m]), 0)
        lr = calib(E[m], gset)
        for k, l in enumerate(LET, 1):
            d = []
            for g in np.nonzero(gset & HIT[:, k - 1])[0]:
                f = (hap[:, g] == k) & (T[:, g] >= MIN_TRUTH) & np.isfinite(lr[:, g])
                b = (hap[:, g] == 0) & (T[:, g] >= MIN_TRUTH) & np.isfinite(lr[:, g])
                if f.sum() >= 1 and b.sum() >= 2:
                    d.append(lr[f, g].mean() - lr[b, g].mean())
            d = np.array(d); med = float(np.median(d))
            want = stored[(stored.cond == 'clean') & (stored.layout == lay) & (stored.method == m) & (stored.founder == l)]
            if abs(round(med, 4) - float(want.median_delta.iloc[0])) > 0.00011 or len(d) != int(want.n_genes.iloc[0]):
                sys.exit(f'reproduction failed {lay} {m} {l}: {med:.4f} vs {want.median_delta.iloc[0]}')
            bs = np.median(d[rng.integers(0, len(d), (2000, len(d)))], axis=1)
            rows.append(dict(layout=lay, method=m, founder=l, median_delta=round(med, 4), ci_lo=round(float(np.percentile(bs, 2.5)), 4),
                             ci_hi=round(float(np.percentile(bs, 97.5)), 4), q25=round(float(np.percentile(d, 25)), 4),
                             q75=round(float(np.percentile(d, 75)), 4), n_genes=len(d)))
    print(lay, 'reproduced', flush=True)

tool = {'hybrid': 'hybrid', 'vf': 'gene probes alone', 'emg': 'EM alone', 'ref': 'ref-unique probes', 'salmon_ref': 'Salmon-ref',
        'salmon_pan': 'Salmon-pan'}
df = pd.DataFrame(rows); df.insert(2, 'tool', df.method.map(tool))
with open(OUT, 'w') as fh:
    fh.write('# file: fig2b_e2_delta_founder_ci.tsv\n'
             '# source: revision/e2_sim/out/est/clean_*.npz, e2_sim/work/design.npz, respacing/probes/probes_s26|77.npz (hit genes), '
             'fast_competitors/e9_prep_cache/ref/genes.json\n'
             '# written by: fig2_recompute_e2_delta.py; derived (recomputed), not copied\n'
             '# note: median_delta and n_genes reproduce e2_sim/out/delta_founder.tsv exactly (checked). ci_lo/ci_hi = gene-bootstrap 95% '
             'percentile interval of the median (2,000 resamples, seed 20261015); q25/q75 = interquartile range of per-gene shifts.\n'
             '# note: founders A C D E F G H = A/J 129S1/SvImJ NOD/ShiLtJ NZO/HlLtJ CAST/EiJ PWK/PhJ WSB/EiJ; clean condition, expected truth.\n')
    df.to_csv(fh, sep='\t', index=False)
print('wrote', OUT, len(df))
