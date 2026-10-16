import json, numpy as np, pandas as pd
REV = '/mnt/project-files/revision'
K, R, SP = 25, 50, 26
G = json.load(open(f'{REV}/fast_competitors/e9_prep_cache/ref/genes.json'))
GNAME = np.array([g['gene_name'] for g in G]); GID = np.array([g['gene_id'] for g in G])
LT = np.array([np.mean([max(L - R + 1, K) for L in g['tx_lens']]) for g in G])
LREF = np.array([g['tx_lens'][0] for g in G])
def meta():
    m = pd.read_csv(f'{REV}/e9_cohort/per_gene.tsv.gz', sep='\t', dtype={'hit': str})
    m['hit'] = m['hit'] == 'True'
    sal = pd.read_csv(f'{REV}/e9_cohort/salmon.tsv.gz', sep='\t').pivot_table(index=['srr', 'gene'], columns='method', values='estimate').reset_index()
    m = m.merge(sal, on=['srr', 'gene'], how='left')
    gi = {}
    for i, n in enumerate(GNAME): gi.setdefault(n, []).append(i)
    m['gi'] = m['gene'].map(lambda n: gi[n][0] if len(gi.get(n, [])) == 1 else -1)
    return m
def gene_matrix(m, est, col):
    """est: (n_samples, n_genes) array in srr order -> column in m."""
    srrs = sorted(m.srr.unique()); si = {s: i for i, s in enumerate(srrs)}
    m[col] = est[m.srr.map(si).to_numpy(), m.gi.to_numpy()]
    return m
