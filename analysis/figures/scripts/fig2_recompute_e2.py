#!/usr/bin/env python3
"""Per-sample error in the read simulation, recomputed for Fig. 2a and Fig. S1.

Reads (read only):  revision/e2_sim/out/est/{cond}_{layout}.npz  (gene x sample estimates, simulation gene order)
                    revision/e2_sim/work/design.npz             (expected counts, haplotypes; same file the scorer used)
Writes:             figures/data/fig2a_e2_per_sample_error.tsv

The metric is the simulation scorer's M2 (revision/e2_sim/scripts/score_e2.py, functions calib and m2), applied per sample instead of
taking the median over the 8 samples: per-sample constant C_s = median(truth / estimate) over B6 sample-genes with truth >= 50,
then the median |log2((C_s * est + 1) / (truth + 1))| over B6 sample-genes with truth >= 50. The median over samples of this value
reproduces M2_err_B6 in revision/e2_sim/out/metrics.tsv; the script checks that and stops if it does not.
Gene sets follow the scorer: core = genes with a finite hybrid value; each method on core genes where it is finite.
"""
import numpy as np, pandas as pd, os, sys

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision/e2_sim'
OUT = '/mnt/project-files/paper/figures/data/fig2a_e2_per_sample_error.tsv'
NS, MIN_TRUTH = 8, 50
D = np.load(f'{REV}/work/design.npz')
NG = len(D['genes'])
B6 = D['hap'].astype(int) == 0
metrics = pd.read_csv(f'{REV}/out/metrics.tsv', sep='\t')


def truth(cond):
    t = np.stack([np.bincount(D['gene_idx'], weights=D['expected'][s], minlength=NG) for s in range(NS)])
    return t * (0.85 if cond == 'premrna' else 1.0)


def per_sample_err(E, T, gset):
    out = []
    for s in range(NS):
        q = B6[s] & gset & (T[s] >= MIN_TRUTH) & (E[s] > 0) & np.isfinite(E[s])
        C = np.median(T[s, q] / E[s, q])
        lr = np.log2((C * E[s] + 1) / (T[s] + 1))
        m = B6[s] & gset & (T[s] >= MIN_TRUTH) & np.isfinite(lr)
        out.append((float(np.median(np.abs(lr[m]))), float(C), int(m.sum())))
    return out


rows = []
for cond, lays in {'clean': ('se50', 'se101', 'pe'), 'premrna': ('pe',), 'gc': ('pe',), 'pos3': ('pe',), 'hex': ('pe',)}.items():
    for lay in lays:
        E = dict(np.load(f'{REV}/out/est/{cond}_{lay}.npz'))
        T = truth(cond)
        core = np.all(np.isfinite(E['hybrid']), 0)
        PK = 'salmon_pan' if 'salmon_pan' in E else 'salmon_pan_flags'
        for m in [k for k in ('hybrid', 'vf', 'emg', 'ref', 'salmon_ref', 'salmon_pan', 'salmon_pan_flags', 'salmon_pan_flags_tpm') if k in E]:
            gset = core & np.all(np.isfinite(E[m]), 0)
            mine = per_sample_err(E[m], T, gset)
            pan = per_sample_err(E[PK], T, gset)          # Salmon-pan on the same genes, as M2_pan_same_genes
            for s in range(NS):
                rows.append(dict(cond=cond, layout=lay, tier='core', method=m, sample=s, err_B6=round(mine[s][0], 6),
                                 C=round(mine[s][1], 4), n_sample_genes=mine[s][2],
                                 pan_err_same_genes=round(pan[s][0], 6), ratio_to_pan=round(mine[s][0] / pan[s][0], 4)))
            med = np.median([x[0] for x in mine])
            ref = metrics[(metrics.cond == cond) & (metrics.layout == lay) & (metrics.method == m) &
                          (metrics.truth == 'expected') & (metrics.tier == 'core')]
            if len(ref):
                want = float(ref.M2_err_B6.iloc[0])
                ok = abs(round(med, 4) - want) <= 0.00011
                print(f'{cond:8s} {lay:6s} {m:17s} median {med:.4f}  metrics.tsv {want:.4f}  {"OK" if ok else "MISMATCH"}')
                if not ok:
                    sys.exit('reproduction check failed')

# extended tier (reported separately; scorer: gset = extended-tier genes with finite v2ext and finite method value)
sg = {g: i for i, g in enumerate(D['genes'].tolist())}
ext = np.zeros(NG, bool)
for line in open('/mnt/project-files/revision/frozen_v2/gene_tiers.tsv'):
    x = line.rstrip('\n').split('\t')
    if x[0] in sg and x[3] == 'extended':
        ext[sg[x[0]]] = True
for lay in ('se50', 'se101', 'pe'):
    E = dict(np.load(f'{REV}/out/est/clean_{lay}.npz')); T = truth('clean')
    for m in ('v2ext', 'salmon_ref', 'salmon_pan'):
        gset = ext & np.all(np.isfinite(E['v2ext']), 0) & np.all(np.isfinite(E[m]), 0)
        mine, pan = per_sample_err(E[m], T, gset), per_sample_err(E['salmon_pan'], T, gset)
        for s in range(NS):
            rows.append(dict(cond='clean', layout=lay, tier='extended', method=m, sample=s, err_B6=round(mine[s][0], 6),
                             C=round(mine[s][1], 4), n_sample_genes=mine[s][2], pan_err_same_genes=round(pan[s][0], 6),
                             ratio_to_pan=round(mine[s][0] / pan[s][0], 4)))
        want = float(metrics[(metrics.cond == 'clean') & (metrics.layout == lay) & (metrics.method == m) & (metrics.truth == 'expected')
                             & (metrics.tier == 'extended')].M2_err_B6.iloc[0])
        med = np.median([x[0] for x in mine]); ok = abs(round(med, 4) - want) <= 0.00011
        print(f'extended {lay:6s} {m:17s} median {med:.4f}  metrics.tsv {want:.4f}  {"OK" if ok else "MISMATCH"}')
        if not ok:
            sys.exit('reproduction check failed')

tool = {'hybrid': 'hybrid', 'v2ext': 'hybrid (extended tier)', 'salmon_pan_flags_tpm': 'Salmon-pan (bias flags, TPM x eff. length, post hoc)', 'vf': 'gene probes alone', 'emg': 'EM alone', 'ref': 'ref-unique probes',
        'salmon_ref': 'Salmon-ref', 'salmon_pan': 'Salmon-pan', 'salmon_pan_flags': 'Salmon-pan (bias flags)'}
df = pd.DataFrame(rows)
df.insert(3, 'tool', df.method.map(tool))
with open(OUT, 'w') as fh:
    fh.write('# file: fig2a_e2_per_sample_error.tsv\n'
             '# source: revision/e2_sim/out/est/*.npz + revision/e2_sim/work/design.npz (expected-count truth)\n'
             '# written by: fig2_recompute_e2.py; derived (recomputed), not copied\n'
             '# note: err_B6 = per-sample median |log2 error| on B6 sample-genes with truth >= 50 after a per-sample constant (simulation metric M2, '
             'per sample). The median over the 8 samples equals M2_err_B6 in e2_sim/out/metrics.tsv for every row (checked by the script).\n'
             '# note: tier=extended rows use the 2,062 extended-tier genes (reported separately). salmon_pan_flags_tpm is Salmon-pan run '
             'with bias flags, counts taken as TPM x effective length.\n'
             '# note: method codes: vf = gene probes alone, emg = EM alone, ref = ref-unique probes, hybrid = PLUCK. pan_err_same_genes = '
             'Salmon-pan on the method\'s genes; ratio_to_pan = err_B6 / pan_err_same_genes per sample.\n')
    df.to_csv(fh, sep='\t', index=False)
print('wrote', OUT, len(df), 'rows')
