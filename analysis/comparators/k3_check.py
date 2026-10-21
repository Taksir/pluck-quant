#!/usr/bin/env python3
"""k3_check.py COHORT SRR OUTDIR_ROOT : pipeline sanity check (not truth).  Spearman, over core genes with stored Salmon-ref >= 50 reads in the same sample,
of kallisto_ref / kallisto_pan / star_fc (>= 0.85) with the stored Salmon-ref gene value; kmerator mean vs kallisto_ref gene TPM (>= 0.70)."""
import os, sys, gzip, json, collections, numpy as np
from scipy.stats import spearmanr
# Folder constants below point to the original working folder; set them for your copy.
coh, srr, root = sys.argv[1:4]; R = '/mnt/project-files/revision'; W = '/home/user/fc'
core = [l.strip() for l in open(f'{W}/km/core_genes.txt') if l.strip()]
def rd(p, cols):
    d = {}
    with gzip.open(p, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); ix = [h.index(c) for c in cols]
        for l in f:
            x = l.rstrip('\n').split('\t'); d[x[0]] = [float(x[i]) for i in ix]
    return d
t2g = dict(l.rstrip('\n').split('\t')[:2] for l in open(f'{W}/fa/tx2gene.tsv'))
sal = collections.defaultdict(float)
if coh == 'keele_kidney':
    with gzip.open(f'{R}/e4_kidney/salmon/{srr}_ref.quant.sf.gz', 'rt') as f:
        next(f)
        for l in f:
            x = l.split('\t'); g = t2g.get(x[0])
            if g: sal[g] += float(x[4])
elif coh == 'keele_liver':
    G = json.load(open(f'{R}/fast_competitors/e9_prep_cache/ref/genes.json')); cnt = collections.Counter(g['gene_name'] for g in G)
    sym2id = {g['gene_name']: g['gene_id'] for g in G if cnt[g['gene_name']] == 1}
    with gzip.open(f'{R}/e9_cohort/salmon.tsv.gz', 'rt') as f:
        next(f)
        for l in f:
            x = l.rstrip('\n').split('\t')
            if x[0] == srr and x[2] == 'salmon_ref' and x[1] in sym2id: sal[sym2id[x[1]]] = float(x[3])
else:
    for l in open(f'{R}/e5_munger/salmon/{srr}_ref.tsv'):
        if l[0] == '#': continue
        x = l.split(); sal[x[0]] = float(x[1])
o = f'{root}/{coh}/{srr}'
kr = rd(f'{o}/kallisto_ref.genes.tsv.gz', ['est_counts', 'tpm']); kp = rd(f'{o}/kallisto_pan.genes.tsv.gz', ['est_counts', 'tpm']) if os.path.exists(f'{o}/kallisto_pan.genes.tsv.gz') else None
sf = rd(f'{o}/featurecounts.genes.tsv.gz', ['count']); km = rd(f'{o}/kmerator.genes.tsv.gz', ['mean'])
sel = [g for g in core if sal.get(g, 0) >= 50]
def sp(a, b, key=lambda v: v[0]): return spearmanr([key(a.get(g, [0.0])) for g in sel], [b[g] for g in sel])[0]
SAL = {g: sal[g] for g in sel}
res = {'cohort': coh, 'srr': srr, 'n_genes': len(sel),
       'kallisto_ref': sp(kr, SAL), 'kallisto_pan': (sp(kp, SAL) if kp is not None else None), 'star_fc': sp(sf, SAL), 'kmerator_vs_kallisto_tpm': sp(km, {g: kr.get(g, [0, 0])[1] for g in sel})}
res['pass'] = res['kallisto_ref'] >= .85 and (res['kallisto_pan'] is None or res['kallisto_pan'] >= .85) and res['star_fc'] >= .85 and res['kmerator_vs_kallisto_tpm'] >= .70
print(json.dumps({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else (bool(v) if isinstance(v, np.bool_) else v)) for k, v in res.items()}))
