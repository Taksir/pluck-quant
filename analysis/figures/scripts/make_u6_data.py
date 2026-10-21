"""Extract the extra source data for Figs S10-S11 and Tables S7, S11, S12.

Reads only revision/ in the original working folder and the existing figure TSVs; writes only
figures/data/ (files prefixed figS10_, tableS7_).
Values are copied from the source files; columns marked 'derived' are computed here.
"""
import json, os, datetime
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision'
OUT = '/mnt/project-files/paper/figures/data'
STAMP = datetime.date(2026, 10, 15).isoformat()


def write(df, name, source, notes=()):
    path = os.path.join(OUT, name)
    with open(path, 'w') as f:
        f.write(f'# file: {name}\n# source: {source}\n')
        f.write(f'# written: {STAMP} by make_u6_data.py; values copied from the source, not recomputed unless marked derived\n')
        for n in notes:
            f.write(f'# note: {n}\n')
        df.to_csv(f, sep='\t', index=False)
    print('wrote', path, len(df))


E6 = f'{REV}/e6_power/out'
DATASETS = [  # key, cis json, label, read length, n
    ('keele_kidney', 'cis_kidney.json', 'Keele kidney', 50),
    ('keele_liver', 'cis_liver.json', 'Keele liver', 50),
    ('e4b_liver', 'cis_e4b_liver.json', 'GSE199702 liver', 101),
    ('e4b_heart', 'cis_e4b_heart.json', 'GSE199702 heart', 101),
    ('e4b_kidney', 'cis_e4b_kidney.json', 'GSE199702 kidney', 101),
]

# --- Fig. S10a: power-model detection curves (mean detections at BH 5% vs number of animals)
rows = []
for key, fn, label, rl in DATASETS:
    d = json.load(open(f'{E6}/{fn}'))
    ns = d['subsample']['n']
    for m, det in d['subsample']['detections'].items():
        for n, v in zip(ns, det):
            rows.append(dict(cohort=label, read_length_bp=rl, method=m, n_animals=n, mean_detections=v))
write(pd.DataFrame(rows), 'figS10_e6_detection_curves.tsv',
      'e6_power/out/cis_{kidney,liver}.json, cis_e4b_{liver,heart,kidney}.json (key subsample)',
      ['method codes: hybrid = PLUCK; salmon_pan, salmon_ref, ref_unique = ref-unique probes; keele_rsem / gbrs = the cohort truth pipeline (not a comparator in the paper; not drawn).'])

# --- Fig. S10b: excess-variance bins and fits D(T) = delta/T + gamma
ex = json.load(open(f'{E6}/excess_fits.json'))
rows = []
labels = {'keele_kidney': ('Keele kidney', 50), 'keele_liver': ('Keele liver', 50), 'munger': ('Munger founders', 100),
          'e4b_liver': ('GSE199702 liver', 101), 'e4b_heart': ('GSE199702 heart', 101), 'e4b_kidney': ('GSE199702 kidney', 101)}
for key, (label, rl) in labels.items():
    e = ex[key]
    for T, D in zip(e['bins_T'], e['bins_D_median']):
        rows.append(dict(cohort=label, read_length_bp=rl, bin_T_median_reads=T, D_median=D,
                         delta_median=e['delta_median'], gamma_median=e['gamma_median']))
write(pd.DataFrame(rows), 'figS10_e6_excess_bins.tsv', 'e6_power/out/excess_fits.json (bins_T, bins_D_median, delta_median, gamma_median)',
      ['D = PLUCK residual variance minus Salmon-pan residual variance (log2^2) in the truth-free scan, median per expression bin; T = gene reads per sample (bin median).',
       'noise_kidney / noise_liver keys (excess against the same-read truth) are not copied: they are the SD-check refit, not the read-factor fit.'])

# --- Fig. S10c and Table S11: read factor for equal noise, tail shares, extra-n
tail_rows = []
for key, (label, rl) in labels.items():
    e = ex[key]
    tail_rows.append(dict(cohort=label, read_length_bp=rl, n_genes=e['n_genes'],
                          pluck_frac_gt2=e['ratio_tail']['frac_gt2'], pluck_frac_gt4=e['ratio_tail']['frac_gt4'],
                          pluck_q50=e['ratio_tail']['q50'], pluck_q95=e['ratio_tail']['q95'],
                          salmon_ref_frac_gt2=e['ratio_tail_salmon_ref']['frac_gt2'],
                          salmon_ref_q50=e['ratio_tail_salmon_ref']['q50']))
write(pd.DataFrame(tail_rows), 'figS10_e6_tail.tsv', 'e6_power/out/excess_fits.json (ratio_tail, ratio_tail_salmon_ref)',
      ['Per-gene residual-variance ratio to Salmon-pan in the truth-free scan; frac_gt2 = share of genes above 2x.'])

model = json.load(open(f'{E6}/model.json'))
rf = [
    # source, layout, delta, gamma, read factor text, read factor low, high, basis
    ('E2 simulation', '50 bp SE', 1.11, 0.0, '1.51', 1.514, 1.514, 'independent truth; model.json depth_factor_sim'),
    ('E2 simulation', '101 bp SE', 0.80, 0.0, '1.38', 1.383, 1.383, 'independent truth; model.json depth_factor_sim'),
    ('E2 simulation', '101 bp PE', 0.57, 0.0, '1.27', 1.274, 1.274, 'independent truth; model.json depth_factor_sim'),
    ('Munger founders', '100 bp SE', ex['munger']['delta_median'], ex['munger']['gamma_median'], '1.38-1.45', 1.38, 1.45, 'truth-free strain ANOVA'),
    ('GSE199702 liver', '101 bp SE', ex['e4b_liver']['delta_median'], ex['e4b_liver']['gamma_median'], '1.57', 1.57, 1.57, 'truth-free cis-eQTL scan'),
    ('GSE199702 heart', '101 bp SE', ex['e4b_heart']['delta_median'], ex['e4b_heart']['gamma_median'], '1.61', 1.61, 1.61, 'truth-free cis-eQTL scan'),
    ('GSE199702 kidney', '101 bp SE', ex['e4b_kidney']['delta_median'], ex['e4b_kidney']['gamma_median'], '1.45', 1.45, 1.45, 'truth-free cis-eQTL scan'),
    ('Keele kidney', '50 bp SE', ex['keele_kidney']['delta_median'], ex['keele_kidney']['gamma_median'], '1.82', 1.82, 1.82, 'truth-free cis-eQTL scan'),
    ('Keele liver', '50 bp SE', ex['keele_liver']['delta_median'], ex['keele_liver']['gamma_median'], '1.8-2.1 (T <= 300)', 1.8, 2.1, 'truth-free cis-eQTL scan; floor gamma > 0 so the factor rises with T'),
]
write(pd.DataFrame(rf, columns=['source', 'layout', 'delta_median', 'gamma_median', 'read_factor_text', 'read_factor_low', 'read_factor_high', 'basis']),
      'figS10_e6_read_factor.tsv', 'e6_power/results.md section 3 table (read factors, text); delta/gamma from e6_power/out/excess_fits.json; simulation delta from results.md section 3',
      ['Read factor = PLUCK needs this many times the reads of Salmon-pan to match its per-gene noise; a multiplier on reads, not an SD ratio.'])

grid = pd.DataFrame(model['extra_n_grid'])
write(grid, 'figS10_e6_extra_n_grid.tsv', 'e6_power/out/model.json (extra_n_grid)',
      ['k = animals PLUCK needs relative to Salmon-pan for one gene with T reads per sample and residual SD resid_sd (log2); extra_n_sim uses the read-simulation alpha, extra_n_real_median_gene the real fit named in real_fit_source.'])

mu = json.load(open(f'{E6}/munger_strain.json'))
write(pd.DataFrame(mu['model_check_by_expression']), 'figS10_e6_munger_model_check.tsv', 'e6_power/out/munger_strain.json (model_check_by_expression)',
      ['Observed vs predicted within-strain variance ratio (PLUCK over Salmon-pan) per expression bin; prediction from the read-simulation alpha.'])

# --- Table S7: every timing block (main-cluster ratios already in fig6c_timing_ratios.tsv; second-cluster runs one per row)
syc = pd.read_csv(f'{REV}/final_comparison/kmer_tools/results/timing/timing_per_run.tsv', sep='\t', comment='#')
write(syc, 'figS11_second_cluster_timing_per_run.tsv', 'final_comparison/kmer_tools/results/timing/timing_per_run.tsv',
      ['Second cluster, AMD EPYC, 4 threads; separate machine from the main cluster; ratios valid only within each paired block. countTags timed with the patched build (sha f39fb662...73e7).'])
