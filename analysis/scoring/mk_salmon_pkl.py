# builds out/salmon.pkl exactly as estimates() does (script's own salmon()), then compares with the stored salmon.pkl.gz
import sys, pickle, gzip, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location('s5', '/mnt/project-files/revision/e5_munger/scripts/score_e5_post_lock_fix.py'); s5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s5)
E = np.load(f'{s5.O}/est.npz', allow_pickle=True); srrs = list(E['srrs'])
sr, sp = s5.salmon(srrs, 'ref'), s5.salmon(srrs, 'pan')
pickle.dump(dict(ref=sr, pan=sp), open(f'{s5.O}/salmon.pkl', 'wb'))
old = pickle.load(gzip.open('/mnt/project-files/revision/e5_munger/salmon.pkl.gz', 'rb'))
ok = all(old[k][s] == v[s] for k, v in (('ref', sr), ('pan', sp)) for s in srrs)
print('salmon.pkl rebuilt for', len(srrs), 'samples; identical to stored salmon.pkl.gz:', ok)
