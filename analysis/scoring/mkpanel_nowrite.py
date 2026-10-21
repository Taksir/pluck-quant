"""Build one deduplicated antisense/both-orientation panel: rc(U) [probe_em EM union incl. s26 gene probes],
rc(probe_broaden sense_code) [v2 extended sets], diagnostic k-mers (dk sites.npz ref/alt) in both orientations.
This copy of mkpanel.py writes only map.npz: the index of every set member in the sorted panel (panel.txt, one 25-mer
per line sorted by 2-bit code, is not written)."""
import numpy as np, sys
# Folder constants below point to the original working folder; set them for your copy.
R = '/mnt/project-files/revision'; W = '/home/user/e4'
K = 25
def rc(c):
    c = c.astype(np.uint64); r = np.zeros_like(c)
    for i in range(K):
        r = (r << np.uint64(2)) | (np.uint64(3) - (c & np.uint64(3))); c = c >> np.uint64(2)
    return r
U = np.load(f'{R}/probe_em/design.npz')['U']
B = np.load(f'{R}/probe_broaden/design.npz')['sense_code']
S = np.load(f'{W}/dk/sites.npz')
dk = np.concatenate([S['ref'], S['alt']]).astype(np.uint64)
aU, aB, d1, d2 = rc(U), rc(B), dk, rc(dk)
pan = np.unique(np.concatenate([aU, aB, d1, d2]))
ix = lambda a: np.searchsorted(pan, a)
np.savez(f'{W}/map.npz', iU=ix(aU), iB=ix(aB), idk_f=ix(d1), idk_r=ix(d2), nsites=len(S['ref']), npanel=len(pan))
print('panel', len(pan), 'U', len(U), 'B', len(B), 'dk sites', len(S['ref']))
