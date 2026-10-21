#!/usr/bin/env python3
"""DO false-cis rate ratio, PLUCK / Salmon-pan, as a plain ratio (no denominator floor).

The scorer (revision/idea1_sim/scripts/score_i1.py) divides by max(Salmon-pan rate, 0.005), which gives 0.488 in DO,
where Salmon-pan's rate (0.33%) is below that floor. This script recomputes the calls from the stored scan p-values
with the scorer's own Benjamini-Hochberg rule, checks the null-call counts against results.json, and reports the
plain ratio with a paired gene bootstrap (1,000 resamples, seed 20261009, as in the scorer).
Reads: revision/idea1_sim/out/design.npz, out/score/scan_p_se50.npz, scan_p_se101.npz, results.json
Writes: figures/data/fig2c_do_fpr_plain.tsv
"""
import json
import numpy as np

# Folder constants below point to the original working folder; set them for your copy.
O = '/mnt/project-files/revision/idea1_sim/out'
OUT = '/mnt/project-files/paper/figures/data/fig2c_do_fpr_plain.tsv'
kind = np.load(f'{O}/design.npz')['kind']
res = json.load(open(f'{O}/score/results.json'))
K0, K1 = kind == 0, kind == 1
fam = K0 | K1


def bh(p, fam):
    idx = np.nonzero(fam & np.isfinite(p))[0]
    o = idx[np.argsort(p[idx])]; m = len(o)
    ok = p[o] <= 0.05 * np.arange(1, m + 1) / m
    c = np.zeros(len(p), bool)
    if ok.any():
        c[o[:np.max(np.nonzero(ok)[0]) + 1]] = True
    return c


rows = ['# file: fig2c_do_fpr_plain.tsv', '# source: idea1_sim/out/score (scan p-values, design kinds); script fig2c_do_fpr_plain.py',
        'layout\tgate\tpluck_null_calls\tpan_null_calls\tn_null\tpluck_rate\tpan_rate\tratio\tci_lo\tci_hi']
for lay in ('se50', 'se101'):
    Z = np.load(f'{O}/score/scan_p_{lay}.npz')
    call, n0m = {}, {}
    for m in ('hybrid', 'salmon_pan'):
        p = Z[f'reads|{m}|do'].astype(float)
        t = fam & np.isfinite(p)
        call[m], n0m[m] = bh(p, t), t & K0
        assert (call[m] & n0m[m]).sum() / n0m[m].sum() == res[lay]['primary'][f'do|{m}']['fpr']
    n0 = n0m['hybrid'] & n0m['salmon_pan']
    a, b = call['hybrid'][n0].astype(float), call['salmon_pan'][n0].astype(float)
    rng = np.random.default_rng(20261009); r = []
    for _ in range(1000):
        i = rng.integers(0, len(a), len(a))
        r.append(a[i].mean() / b[i].mean())
    lo, hi = np.percentile(r, [2.5, 97.5])
    rows.append(f'{lay}\tG1a_fpr_do\t{int(a.sum())}\t{int(b.sum())}\t{len(a)}\t{a.mean():.6f}\t{b.mean():.6f}\t{a.mean() / b.mean():.3f}\t{lo:.3f}\t{hi:.3f}')
open(OUT, 'w').write('\n'.join(rows) + '\n')
print('\n'.join(rows))
