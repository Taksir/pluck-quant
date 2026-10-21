# Recompute the Munger core estimates (hybrid, vf gene probes, ref-unique, emg, capture) from the stored counts with the scoring script's own
# functions (seg, gene_probe, em_model: unchanged) and compare with the stored est.npz.  The extended tier is not recomputed.
import gzip, importlib.util, io, numpy as np
spec = importlib.util.spec_from_file_location('s5', '/mnt/project-files/revision/e5_munger/scripts/score_e5_post_lock_fix.py'); s5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s5)
E = np.load('/mnt/project-files/revision/e5_munger/est.npz', allow_pickle=True); srrs = [str(s) for s in E['srrs']]
s5.counts = lambda s: np.load(io.BytesIO(gzip.open(f'/mnt/project-files/revision/e5_munger/counts/{s}.npy.gz').read()))
C = [s5.counts(s) for s in srrs]
Uc = np.stack([s5.seg(c, 'U76') for c in C]); Pc = np.stack([s5.seg(c, 's76') for c in C])
S_vf, vf, m_vf = s5.gene_probe(Pc, 'vf'); _, ref, m_ref = s5.gene_probe(Pc, 'ref')
em, M, tg_i, tg_u, gmap, n_emg = s5.em_model(); em.R = s5.R
emg = np.full((len(srrs), s5.NGJ), np.nan)
for i, c in enumerate(Uc):
    X, it = M.fit(c, bias=False); gsum = np.bincount(tg_i, X, len(tg_u)); ok = gmap >= 0; emg[i, gmap[ok]] = gsum[ok]
fo = np.nansum(S_vf, 0) / np.maximum(np.nansum(emg, 0), 1e-9)
hyb = np.where(fo[None, :] >= s5.TAU, vf, emg)
for k, v in dict(hybrid=hyb, vf=vf, ref=ref, emg=emg, capture=fo, m_vf=m_vf, m_ref=m_ref, n_emg=n_emg).items():
    o = np.asarray(E[k], float); n = np.asarray(v, float)
    same_nan = np.array_equal(np.isnan(o), np.isnan(n)); d = np.nanmax(np.abs(o - n)) if o.size else 0
    print(k, 'nan-pattern equal:', same_nan, 'max abs diff:', d, 'max rel diff:', np.nanmax(np.abs(o - n) / np.maximum(np.abs(o), 1e-9)))
