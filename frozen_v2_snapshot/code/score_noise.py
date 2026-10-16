# Per-gene residual noise across strains: lr = log2((C est + 1)/(truth + 1)), C per sample on B6 genes (as score_e9).
# Robust within-gene SD = 1.4826 * MAD of lr across samples (genes with >= 10 B6/nonB6 samples with truth >= 50).
import pandas as pd, numpy as np, json
sal = pd.read_csv('score/salmon.tsv', sep='\t').pivot_table(index=['srr', 'gene'], columns='method', values='estimate').reset_index()
out = {}
for sp in (26, 51, 77, 101, 126):
    w = pd.read_csv(f'score/s{sp}/per_gene.tsv', sep='\t').merge(sal, on=['srr', 'gene'])
    w = w[w.stratum.isin(['B6', 'nonB6'])]
    meths = {'ref': 'est_ref', 'vf': 'est_vf', 'salmon_ref': 'salmon_ref', 'salmon_pan': 'salmon_pan'}
    for m, col in meths.items():
        b = w[(w.stratum == 'B6') & (w.truth >= 50) & (w[col] > 0)]
        C = (b.truth / b[col]).groupby(b.srr).median()
        w['lr_' + m] = np.log2((w.srr.map(C) * w[col] + 1) / (w.truth + 1))
    x = w[w.truth >= 50]
    g = x.groupby('gene').filter(lambda d: len(d) >= 10).groupby('gene')
    res = {}
    for m in meths:
        mad = g['lr_' + m].agg(lambda v: 1.4826 * np.median(np.abs(v - np.median(v))))
        res[m] = dict(median_within_gene_sd=round(float(mad.median()), 4), median_abs_err=round(float(np.median(np.abs(x['lr_' + m]))), 4), n_genes=int(len(mad)))
    out[sp] = res
    print(sp, json.dumps(res))
json.dump(out, open('score/noise.json', 'w'), indent=1)
