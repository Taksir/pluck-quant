"""Design compile for the hybrid estimator: turn design.npz (+genome_count, gc, gene-probe set) into the static arrays that the
per-sample step needs, so a sample run needs only numpy.  Same EM set rule as run_emg.py (emg) and dump.py.
usage: compile_design.py DESIGN_DIR OUT_DIR R PROBES_NPZ
  DESIGN_DIR has design.npz genome_count.npy gc.npy ; R = read length (50 or 101)"""
import sys, os, json, time, numpy as np
# Folder constants below point to the original working folder; set them for your copy.
sys.path.insert(0, os.environ.get('FC_CODE', '/home/user/e3/code'))
ddir, out, R, probes = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]   # probes = gene-probe set (probes_s26.npz at 50 bp, probes_s77.npz at 101 bp)
os.makedirs(out, exist_ok=True); t0 = time.time()
import em
em.R = R
K = 25
D = dict(np.load(f'{ddir}/design.npz')); gc = np.load(f'{ddir}/gc.npy'); gcount = np.load(f'{ddir}/genome_count.npy')
tg_u, tg_i = np.unique(D['tgene'], return_inverse=True); tg_i = tg_i.ravel()
pairs = np.unique(np.stack([D['cp'], tg_i[D['ct']]], 1), axis=0); ngp = np.bincount(pairs[:, 0], minlength=len(D['U']))
use = D['inB'] & (gcount <= np.maximum(ngp, 1)); use_t = use & (ngp == 1)
M = em.Model(D, use_t, gc)
h5 = np.ones(em.NB); h3 = np.ones(em.NB); hg = np.ones(8)
a = M.rates(h5, h3, hg); A = np.bincount(M.kt, weights=a, minlength=M.ntx)
m = M.ec_of_probe >= 0
nprb = np.bincount(M.ec_of_probe[m], minlength=M.nec).astype(np.float64)
ea = a[M.e_ix] * nprb[M.e_ec]; et = M.kt[M.e_ix].astype(np.int32)
sys.path.insert(0, os.environ.get('FC_CODE', '/home/user/e3/code')); import e9rs
Pg = np.load(probes); vfm = Pg['in_vf']
Gcodes = np.array([e9rs.seqcodes(e9rs.enc(k))[0][0] for k in Pg['kmer'][vfm].tolist()], np.uint64); ag = Pg['gene'][vfm].astype(np.int64)
G = json.load(open(os.environ.get('FC_GENES', '/mnt/project-files/revision/fast_competitors/e9_prep_cache/ref/genes.json')))
GID = np.array([g['gene_id'] for g in G]); ng = len(G)
LT = np.array([np.mean([max(L - R + 1, K) for L in g['tx_lens']]) for g in G])
mg = np.bincount(ag, minlength=ng)
gpos = {g: i for i, g in enumerate(GID)}; gmap = np.array([gpos.get(g, -1) for g in tg_u])
em_idx = np.nonzero(m)[0]
def dec(c):
    c = c.astype(np.uint64); res = np.zeros((len(c), K), np.uint8); x = c.copy()
    for i in range(K - 1, -1, -1):
        res[:, i] = x & np.uint64(3); x >>= np.uint64(2)
    return res
U = D['U']
em_codes = U[em_idx]
pcodes = np.unique(np.concatenate([em_codes, Gcodes]))   # sorted; equals the old sel panel when the gene probes are inside U
b = dec(pcodes); comp = (3 - b)[:, ::-1]   # antisense = reverse complement
chars = np.frombuffer(b'ACGT', np.uint8)[comp]
with open(f'{out}/panel.txt', 'wb') as f:
    f.write(b'\n'.join(bytes(r) for r in chars) + b'\n')
pos = np.searchsorted(U, pcodes); pos[pos == len(U)] = 0; inU = U[pos] == pcodes
uidx = pos[inU]; uline = np.nonzero(inU)[0]
gline = np.searchsorted(pcodes, Gcodes)
np.savez(f'{out}/static.npz', ec_of_probe=M.ec_of_probe.astype(np.int32), e_ec=M.e_ec.astype(np.int32), et=et, ea=ea, A=A,
         gof=tg_i.astype(np.int64), gline=gline, ag=ag, LT=LT, mg=mg, gmap=gmap, uidx=uidx, uline=uline, pcodes=pcodes, npanel=np.int64(len(pcodes)),
         nU=np.int64(len(U)), ng=np.int64(ng), nec=np.int64(M.nec), ntx=np.int64(M.ntx), R=np.int64(R), ntg=np.int64(len(tg_u)))
np.save(f'{out}/GID.npy', GID)
print('compiled', round(time.time() - t0, 1), 's; union', len(U), 'panel', len(pcodes), 'gene probes', len(Gcodes), 'EM probes', len(em_idx), 'ECs', M.nec, 'entries', len(et), flush=True)
