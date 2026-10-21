#!/usr/bin/env python3
"""Table S8: metrics of scoring set A (all methods) for all, unflagged and flagged genes under both per-sample calibrations,
for every pipeline of the set. 'as scored' = constants fitted inside the subset being scored (as in the scorer);
'all-gene' = constants of the all-gene score. Both constant tables are computed before the calibration function is replaced.
Usage: python3 tableS8_calib_u5.py <cohort> <out.tsv>   (needs scipy; a few minutes per cohort)."""
import importlib.util, sys, json, os
# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision/final_comparison'; coh, out = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location('fc', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scoring', 'fc_gather_score.py')); fc = importlib.util.module_from_spec(spec); spec.loader.exec_module(fc)
fc.load_scorer(); c = fc.load_cohort(coh, f'{REV}/frames')
M = fc.EXISTING + fc.acc_for(coh)                       # method list of set A, as in fc_gather_score.py cmd_score
fc.attach_tools(c, f'{REV}/results_ll/acc', fc.acc_for(coh))
fc._CTX[coh] = c; sets, _ = fc.flag_sets(c)
fc.set_calibration(coh); cal = fc.S9.calibrate
C = {'all-gene': cal(c.recs, M)}
own = {lab: cal({k: r for k, r in c.recs.items() if k[1] in sets[lab]}, M) for lab in ('flagged', 'unflagged')}
fc.set_calibration = lambda cohort: None
cols = ['n_delta_genes', 'median_delta', 'delta_ci95', 'spurious_rate', 'within_gene_spearman', 'n_genes_spearman', 'noise_sd', 'n_genes_noise']
with open(out, 'w') as fh:
    fh.write('cohort\tsubset\tcalibration\tmethod\t' + '\t'.join(cols) + '\n')
    jobs = [('all', 'as scored', C['all-gene'])] + [(lab, cl, k) for lab in ('unflagged', 'flagged') for cl, k in (('as scored', own[lab]), ('all-gene', C['all-gene']))]
    for lab, cl, K in jobs:
        fc.S9.calibrate = lambda recs, methods, K=K: K
        res = fc.score_frame((coh, lab, 'A', M, sets[lab]))
        for m in M:
            x = res['metrics'][m]; x['n_delta_genes'] = res['n_delta_genes']
            fh.write(f'{coh}\t{lab}\t{cl}\t{m}\t' + '\t'.join(str(x.get(k)) for k in cols) + '\n')
        fh.flush(); print(coh, lab, cl, 'done', flush=True)
