#!/usr/bin/env python3
"""Copy the planted cis-eQTL simulation values that Fig. 2c and Fig. S3 need and that fig2c_idea1_gates.tsv does not carry
(per-pipeline false-cis rates and detections in CC and DO, per-founder shifts of every pipeline, effect-size shifts,
the null-relative check, the recalibration on all null core genes and the branch split). Values are copied, not recomputed.

Reads (read only): revision/idea1_sim/out/score/{results,null_relative_bias,posthoc_calib_results,posthoc_g5}.json
Writes: figures/data/fig2c_idea1_rates.tsv, fig2c_idea1_power.tsv, figS3_idea1_founder.tsv, figS3_idea1_effect.tsv,
figS3_idea1_branch.tsv, figS3_idea1_pairs.tsv
"""
import json
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
S = '/mnt/project-files/revision/idea1_sim/out/score'
O = '/mnt/project-files/paper/figures/data'
R = json.load(open(f'{S}/results.json'))
NR = json.load(open(f'{S}/null_relative_bias.json'))
PC = json.load(open(f'{S}/posthoc_calib_results.json'))
PG = json.load(open(f'{S}/posthoc_g5.json'))
TOOL = {'hybrid': 'hybrid', 'vf': 'gene probes alone', 'emg': 'EM alone', 'ref': 'ref-unique probes', 'salmon_ref': 'Salmon-ref',
        'salmon_pan': 'Salmon-pan', 'truth': 'truth (expected reads)'}
HEAD = ('# file: {f}\n# source: revision/idea1_sim/out/score/{src}\n'
        '# written by: fig2_extract_idea1.py; values copied from the source, not recomputed\n'
        '# note: method codes: hybrid = PLUCK, vf = gene probes alone, emg = EM alone, ref = ref-unique probes. One seed (20261008).\n')


def write(df, f, src, note=''):
    with open(f'{O}/{f}', 'w') as fh:
        fh.write(HEAD.format(f=f, src=src) + (f'# note: {note}\n' if note else ''))
        df.to_csv(fh, sep='\t', index=False)
    print('wrote', f, len(df))


rows = []
for lay in ('se50', 'se101'):
    for key, v in R[lay]['primary'].items():
        pop, m = key.split('|')
        rows.append(dict(population=pop, layout=lay, method=m, tool=TOOL[m], n_tested=v['n_tested'], n_null=v['n_null'],
                         n_eqtl=v['n_eqtl'], calls=v['calls'], null_calls=v['null_calls'], eqtl_calls=v['eqtl_calls'], fpr=v['fpr'],
                         fdp=v['fdp'], power=v['power'], lambda_gc_null=v['lambda_gc_null'],
                         **{f'power_b{b}': v[f'power_b{b}'] for b in ('0.25', '0.5', '0.75', '1')}))
rates = pd.DataFrame(rows)
# detections ratios: G2 (CC, compared with a threshold) and the reported-only DO ratio, copied as stored (they count genes tested by both methods,
# so they differ slightly from eqtl_calls in the primary rows)
pw = []
for lay in ('se50', 'se101'):
    for pop, st in (('cc', R[lay]['gates']['G2_power_cc']), ('do', R[lay]['reported_only']['G2_power_do'])):
        pw.append(dict(population=pop, layout=lay, hybrid_detected=st['hybrid_detected'], pan_detected=st['pan_detected'],
                       ratio=st['ratio'], threshold=st['threshold'], gated='yes' if pop == 'cc' else 'no (reported only)',
                       **{f'power_b{b}_{w}': st['per_level'][b][w] for b in ('0.25', '0.5', '0.75', '1') for w in ('hybrid', 'pan')}))
write(pd.DataFrame(pw), 'fig2c_idea1_power.tsv', 'results.json (gates.G2_power_cc, reported_only.G2_power_do)',
      'detections = planted cis-eQTL genes called at BH 5%; ratio = PLUCK / Salmon-pan (limit 0.90, CC only).')
write(rates, 'fig2c_idea1_rates.tsv', 'results.json (se50|se101 -> primary)',
      'fpr = null genes called at BH 5% / null genes tested (false cis-eQTL rate). For detections ratios use fig2c_idea1_power.tsv.')

rows = []
for lay in ('se50', 'se101'):
    for m, v in R[lay]['bias'].items():
        for f, x in v['founder'].items():
            rows.append(dict(layout=lay, method=m, tool=TOOL[m], founder=f, median_delta=x['median_delta'], n=x['n'],
                             frac_beyond_bound=x['frac_beyond_bound'], q90_abs_delta=x['q90_abs_delta'], source='locked'))
    for m, v in NR[lay].items():
        for f, x in v.items():
            rows.append(dict(layout=lay, method=m, tool=TOOL[m], founder=f, n=x['n'], frac_beyond_bound=x['frac_beyond_bound'],
                             frac_beyond_bound_B6split_null=x['frac_beyond_bound_B6split_null'], excess=x['excess'],
                             excess_se=x['excess_se'], source='null-relative check (reported only)'))
    for m, v in PC[lay]['bias'].items():
        for f, x in v['founder'].items():
            rows.append(dict(layout=lay, method=m, tool=TOOL[m], founder=f, median_delta=x['median_delta'], n=x['n'],
                             frac_beyond_bound=x['frac_beyond_bound'], q90_abs_delta=x['q90_abs_delta'],
                             source='post hoc recalibration on all null core genes'))
write(pd.DataFrame(rows), 'figS3_idea1_founder.tsv', 'results.json (bias), null_relative_bias.json, posthoc_calib_results.json (bias)',
      'CC set, hit genes per founder; frac_beyond_bound = share of genes beyond 0.05 log2 by the noise-aware test (check G4b, limit 0.10).')

rows = []
for lay in ('se50', 'se101'):
    for src, D in (('locked', R), ('post hoc recalibration', PC)):
        for m, v in D[lay]['effect_shift'].items():
            for b, x in v.items():
                rows.append(dict(layout=lay, method=m, tool=TOOL[m], beta=float(b), median_shift=x['median_shift'], n=x['n'], source=src))
        for m, v in D[lay]['bias'].items():
            rows.append(dict(layout=lay, method=m, tool=TOOL[m], beta=None, median_shift=None, n=None, source=src,
                             M2_null_cc=v.get('M2_null_cc')))
write(pd.DataFrame(rows), 'figS3_idea1_effect.tsv', 'results.json and posthoc_calib_results.json (effect_shift, bias.M2_null_cc)',
      'median_shift = effect-size bias at each planted |beta| (check G3, CC, limit 0.05); M2_null_cc = context statistic (median |log2 error| on '
      'null genes), for each calibration in the source column.')

rows = []
for lay in ('se50', 'se101'):
    for sw in ('swM', 'swR'):
        for k, v in PG[lay][sw].items():
            if isinstance(v, dict):
                branch, m = k.split('|')
                rows.append(dict(layout=lay, switch=sw, branch=branch, method=m, tool=TOOL[m], n=v['n'], call_rate=v['call_rate']))
        rows.append(dict(layout=lay, switch=sw, branch='share of switch genes on gene-probe branch', method='hybrid', tool='hybrid',
                         call_rate=PG[lay][sw]['share_genes_on_gene_probe_branch']))
write(pd.DataFrame(rows), 'figS3_idea1_branch.tsv', 'posthoc_g5.json',
      'switch-gene call rate (p < 1e-3, CC) split by the branch PLUCK used for the gene; Salmon-pan values in the matched unit.')

rows = []
for lay in ('se50', 'se101'):
    for src, D in (('locked', R), ('post hoc recalibration', PC)):
        g = D[lay]['gates']['G4c_pair_median_do']
        for pair, v in g['values'].items():
            rows.append(dict(layout=lay, pair=pair, median_delta=v, n=g['n'][pair], source=src))
write(pd.DataFrame(rows), 'figS3_idea1_pairs.tsv', 'results.json and posthoc_calib_results.json (gates.G4c_pair_median_do)',
      'DO set: median shift of PLUCK for each of the 28 heterozygous founder pairs (check G4c, limit 0.02).')
