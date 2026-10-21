# EM variant: drop probes lying within 100 nt of either end of any compatible transcript (same END_MARGIN as the gene probes)
import sys, os, time, numpy as np
from common import *
import em
D = dict(np.load('design.npz')); gc = np.load('gc.npy'); gcount = np.load('genome_count.npy')
tg_u, tg_i = np.unique(D['tgene'], return_inverse=True); tg_i = tg_i.ravel(); em.GENE_OF_TX = tg_i
pairs = np.unique(np.stack([D['cp'], tg_i[D['ct']]], 1), axis=0); ngp = np.bincount(pairs[:, 0], minlength=len(D['U']))
use = D['inB'] & (gcount <= np.maximum(ngp, 1))
L = D['tlen'][D['ct']]; nearend = (D['cpos'] < 100) | (D['cpos'] > L - 100 - K)
bad = np.zeros(len(use), bool); bad[D['cp'][nearend]] = True
use_t = use & (ngp == 1)
print('EM probes gene-unique', int(use_t.sum()), 'of', int(use.sum()), flush=True)
M = em.Model(D, use_t, gc)
gpos = {g: i for i, g in enumerate(GID)}; gmap = np.array([gpos.get(g, -1) for g in tg_u])
for s in sys.argv[1:]:
    f = f'emest/{s}_emg.npz'
    if os.path.exists(f): continue
    c = np.load(f'counts/{s}.npy').astype(np.float64)
    t0 = time.time(); X, it = M.fit(c, bias=False); dt = time.time() - t0
    gsum = np.bincount(tg_i, X, len(tg_u)); est = np.full(len(GID), np.nan); ok = gmap >= 0; est[gmap[ok]] = gsum[ok]
    np.savez_compressed(f, gene=est, X=X.astype(np.float32), sec=dt, it=it)
    print(s, 'emg', round(dt, 1), it, flush=True)
