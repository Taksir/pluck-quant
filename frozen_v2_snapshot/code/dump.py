# Dump EC arrays (uniform model) for one sample for the C EM; compare with the numba result.
import sys, numpy as np
sys.path.insert(0, '/home/user/pe'); import os; os.chdir('/home/user/pe')
from common import *
import em
D = dict(np.load('design.npz')); gc = np.load('gc.npy'); gcount = np.load('genome_count.npy')
tg_u, tg_i = np.unique(D['tgene'], return_inverse=True); tg_i = tg_i.ravel()
pairs = np.unique(np.stack([D['cp'], tg_i[D['ct']]], 1), axis=0); ngp = np.bincount(pairs[:, 0], minlength=len(D['U']))
use = D['inB'] & (gcount <= np.maximum(ngp, 1))
if len(sys.argv) > 2 and sys.argv[2] == 'emg': use = use & (ngp == 1)
o = 'cem/in_' + (sys.argv[2] if len(sys.argv) > 2 else 'emu')
M = em.Model(D, use, gc)
s = sys.argv[1]; c = np.load(f'counts/{s}.npy').astype(np.float64)
h5 = np.ones(em.NB); h3 = np.ones(em.NB); hg = np.ones(8)
a = M.rates(h5, h3, hg); A = np.bincount(M.kt, weights=a, minlength=M.ntx)
m = M.ec_of_probe >= 0
cec = np.bincount(M.ec_of_probe[m], weights=c[m], minlength=M.nec); nprb = np.bincount(M.ec_of_probe[m], minlength=M.nec).astype(float)
ea = a[M.e_ix] * nprb[M.e_ec]; et = M.kt[M.e_ix].astype(np.int64)
keep = cec[M.e_ec] > 0  # zero-count classes add nothing to the E-step; their rates are already in A
nz = np.nonzero(cec > 0)[0]; remap = -np.ones(M.nec, np.int64); remap[nz] = np.arange(len(nz))
e_ec = remap[M.e_ec[keep]]; ea = ea[keep]; et = et[keep]; cec = cec[nz]
eoff = np.searchsorted(e_ec, np.arange(len(nz) + 1)).astype(np.int64)
os.makedirs(o, exist_ok=True)
for n, x in (('eoff', eoff), ('et', et), ('ea', ea.astype(np.float64)), ('cec', cec), ('A', A), ('gof', tg_i.astype(np.int64))):
    x.tofile(f'{o}/{n}.bin')
print('nonzero ECs', int((cec > 0).sum()), 'of', M.nec)
