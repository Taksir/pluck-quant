"""Per-sample step after the counter: counts -> equivalence-class arrays -> C EM -> gene estimates -> hybrid switch.
usage: run_sample.py STATIC_DIR COUNTS.txt OUT_PREFIX [tau]
Writes OUT_PREFIX.npy (estimate per gene, genes.json order), prints stage times (stderr)."""
import sys, os, time, subprocess, numpy as np
t0 = time.time()
sd, cf, outp = sys.argv[1:4]; tau = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
S = np.load(f'{sd}/static.npz')
t1 = time.time()
import pandas as pd; c_sel = pd.read_csv(cf, header=None, dtype=np.int64, engine='c').iloc[:, 0].to_numpy().astype(np.float64)
c_sel = c_sel[:int(S['npanel'])]   # v2 panels append extra probes after the v1 panel
c = np.zeros(int(S['nU'])); c[S['uidx']] = c_sel[S['uline']]
ec = S['ec_of_probe']; m = ec >= 0
nec = int(S['nec'])
cec = np.bincount(ec[m], weights=c[m], minlength=nec)
e_ec = S['e_ec']; keep = cec[e_ec] > 0
nz = np.nonzero(cec > 0)[0]; remap = -np.ones(nec, np.int64); remap[nz] = np.arange(len(nz))
e_ec2 = remap[e_ec[keep]]; ea = S['ea'][keep]; et = S['et'][keep].astype(np.int64); cec2 = cec[nz]
eoff = np.searchsorted(e_ec2, np.arange(len(nz) + 1)).astype(np.int64)
wd = f'{outp}_em'; os.makedirs(wd, exist_ok=True)
for n, x in (('eoff', eoff), ('et', et), ('ea', ea.astype(np.float64)), ('cec', cec2), ('A', S['A']), ('gof', S['gof'])):
    x.tofile(f'{wd}/{n}.bin')
t2 = time.time()
r = subprocess.run([os.environ.get('FC_EM', '/home/user/e3/bin2/em'), wd, '1e-3'], capture_output=True, text=True)
t3 = time.time()
X = np.fromfile(f'{wd}/X.bin', dtype=np.float64)
ntg = int(S['ntg']); ng = int(S['ng'])
gsum = np.bincount(S['gof'], X, ntg); gmap = S['gmap']
em_est = np.full(ng, np.nan); ok = gmap >= 0; em_est[gmap[ok]] = gsum[ok]
Sg = np.bincount(S['ag'], weights=c_sel[S['gline']], minlength=ng); mg = S['mg']
Ng = np.where(mg > 0, Sg * S['LT'] / (np.maximum(mg, 1) * (int(S['R']) - 25 + 1)), np.nan)
f = np.where(em_est > 0, Ng / em_est, np.nan)       # single-sample observed capture (cohort sum in a cohort run)
est = np.where((mg > 0) & (f >= tau), Ng, em_est)
est = np.where(np.isnan(est), em_est, est)
np.save(f'{outp}.npy', est)
t4 = time.time()
import shutil; shutil.rmtree(wd, ignore_errors=True)
print(f'load_static {t1-t0:.2f} counts_to_ec {t2-t1:.2f} em {t3-t2:.2f} hybrid {t4-t3:.2f} total {t4-t0:.2f} | {r.stdout.strip()[:120]}', file=sys.stderr)
