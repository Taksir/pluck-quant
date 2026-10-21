"""Build Tables S7, S11 and S12 (Markdown on stdout, TSV for the supplement).

Reads only figures/data; prints Markdown to stdout and writes tableS7_timing_all_blocks.tsv.
"""
import sys
# Folder constants below point to the original working folder; set them for your copy.
sys.path.insert(0, '/mnt/project-files/paper/figures/scripts')
import pandas as pd
from u6_common import read, CODE, DATA

LAB = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags',
       'Salmon-pan 25%', 'Salmon-pan 10%', 'Salmon-ref 25%', 'Salmon-ref 10%']
SMP = [('50', 'CC021', 'CC021, 50 bp'), ('50', 'CC001', 'CC001, 50 bp'), ('101', 'LiverCC001F', 'Liver CC001 F, 101 bp')]

tr = read('fig6c_timing_ratios.tsv')
tr['label'] = tr.method.map(CODE)
reg = tr[tr.kind == 'regular']

out = []
P = out.append
P('### Table S7. Wall time, CPU time and memory per sample in every timed block\n')
for T in (1, 4, 8):
    P(f'**{T} thread{"s" if T > 1 else ""}' + (' (one run per pipeline)' if T == 1 else ' (median of 3 repetitions; min to max in brackets)') + '**\n')
    hdr = '| Pipeline | ' + ' | '.join(f'{s[2]}: wall (s) | ratio to PLUCK' for s in SMP) + ' |'
    P(hdr)
    P('|' + '---|' * (1 + 2 * len(SMP)))
    for lab in LAB:
        cells = []
        for lay, smp, _ in SMP:
            r = reg[(reg.layout.astype(str) == lay) & (reg['sample'] == smp) & (reg['T'] == T) & (reg.label == lab)]
            if r.empty:
                cells += ['not timed', '']
                continue
            r = r.iloc[0]
            cells.append(f'{r.wall_median_s:.1f}' if T == 1 else f'{r.wall_median_s:.1f} [{r.wall_min:.1f}, {r.wall_max:.1f}]')
            if lab == 'PLUCK':
                cells.append('1')
            elif T == 1:
                cells.append(f'{r.ratio_median:.1f}')
            else:
                cells.append(f'{r.ratio_median:.1f} [{r.ratio_min:.1f}, {r.ratio_max:.1f}]')
        P(f'| {lab} | ' + ' | '.join(cells) + ' |')
    P('')

P('**CPU-time ratio to PLUCK and peak memory, 4 threads (medians of 3 repetitions)**\n')
P('| Pipeline | ' + ' | '.join(f'{s[2]}: CPU ratio | peak RSS (MB)' for s in SMP) + ' |')
P('|' + '---|' * (1 + 2 * len(SMP)))
for lab in LAB:
    cells = []
    for lay, smp, _ in SMP:
        r = reg[(reg.layout.astype(str) == lay) & (reg['sample'] == smp) & (reg['T'] == 4) & (reg.label == lab)]
        if r.empty:
            cells += ['not timed', '']
            continue
        r = r.iloc[0]
        cells += ['1' if lab == 'PLUCK' else f'{r.cpu_ratio_median:.1f}', f'{r.rss_mb_median:,.0f}']
    P(f'| {lab} | ' + ' | '.join(cells) + ' |')
P('')

P('**CPU-time ratio to PLUCK at 1 and 8 threads (1 thread: one run; 8 threads: median of 3 repetitions)**\n')
P('| Pipeline | ' + ' | '.join(f'{s[2]}: 1 thread | 8 threads' for s in SMP) + ' |')
P('|' + '---|' * (1 + 2 * len(SMP)))
for lab in LAB[1:]:
    cells = []
    for lay, smp, _ in SMP:
        for T in (1, 8):
            r = reg[(reg.layout.astype(str) == lay) & (reg['sample'] == smp) & (reg['T'] == T) & (reg.label == lab)]
            cells.append('not timed' if r.empty else f'{r.iloc[0].cpu_ratio_median:.1f}')
    P(f'| {lab} | ' + ' | '.join(cells) + ' |')
P('')

gb = tr[tr.kind == 'gbrs']
P('**GBRS blocks (one repetition each; timing only, accuracy not compared)**\n')
P('| Sample, threads | GBRS wall (s) | PLUCK wall (s) | GBRS / PLUCK | Salmon-pan / PLUCK | Salmon-ref / PLUCK | GBRS peak RSS (MB) |')
P('|---|---|---|---|---|---|---|')
for (lay, smp, T), g in gb.groupby(['layout', 'sample', 'T'], sort=False):
    gg = {r.label: r for r in g.itertuples()}
    name = dict((s[1], s[2]) for s in SMP)[smp]
    P(f'| {name}, {T} | {gg["GBRS"].wall_median_s:,.1f} | {gg["PLUCK"].wall_median_s:.2f} | {gg["GBRS"].ratio_median:.1f} | '
      f'{gg["Salmon-pan"].ratio_median:.1f} | {gg["Salmon-ref"].ratio_median:.1f} | {gg["GBRS"].rss_mb_median:,.0f} |')
P('')

sy = read('fig6c_rnaskim_timing_second_cluster.tsv')
P('**Second-cluster blocks (AMD EPYC, 4 threads; a different machine from the main-cluster node, so ratios hold only within each block)**\n')
P('| Sample | RNA-Skim wall (s) | Compared pipeline | its wall (s) | RNA-Skim / pipeline, median [min, max] |')
P('|---|---|---|---|---|')
for r in sy.itertuples():
    name = dict((s[1], s[2]) for s in SMP)[r.sample]
    P(f'| {name} | {r.median_wall_rnaskim_s:.1f} | {CODE[r.other]} | {r.median_wall_other_s:.2f} | {r.median_ratio_rnaskim_over_other:.2f} [{r.min_ratio:.2f}, {r.max_ratio:.2f}] |')
P('')
rss = read('fig6c_rnaskim_peak_rss.tsv')
pv = rss.pivot(index='method', columns='sample', values='rss_mb')
P('Peak RSS in the second-cluster blocks (MB): ' + '; '.join(
    f'{CODE[m]} ' + ' / '.join(f'{pv.loc[m, s]:,.0f}' for s in ['CC021', 'CC001', 'LiverCC001F'])
    for m in ['hybrid', 'salmon_ref', 'kmerator_counttags', 'rnaskim']) + ' (CC021 / CC001 / Liver CC001 F).\n')

e3 = read('fig6c_e3_cloud_paired_speedups.tsv')
P('**Cloud VM replication (4 vCPU, 15 GB; paired back-to-back runs; core tier = v1, core + coverage tier = v2)**\n')
P('| Sample | Paired rep | Salmon-ref (s) | Salmon-pan (s) | core (s) | core + coverage (s) | Salmon-ref / core | Salmon-ref / core + coverage | Salmon-pan / core | Salmon-pan / core + coverage |')
P('|---|---|---|---|---|---|---|---|---|---|')
f = lambda v, d=2: '' if pd.isna(v) else f'{v:.{d}f}'
for r in e3.itertuples():
    P(f'| {r.sample}, {r.layout} bp | {r.rep} | {f(r.salmon_ref)} | {f(r.salmon_pan)} | {f(r.v1)} | {f(r.v2)} | {f(r.ref_v1)} | {f(r.ref_v2)} | {f(r.pan_v1)} | {f(r.pan_v2)} |')
P('')

# machine-readable S7: all main-cluster rows (regular + GBRS blocks)
keep = tr[['kind', 'layout', 'sample', 'T', 'label', 'n_reps', 'wall_median_s', 'wall_min', 'wall_max', 'ratio_median',
           'ratio_min', 'ratio_max', 'cpu_ratio_median', 'rss_mb_median']].rename(columns={'label': 'pipeline'})
with open(f'{DATA}/tableS7_timing_all_blocks.tsv', 'w') as fh:
    fh.write('# file: tableS7_timing_all_blocks.tsv\n# source: fig6c_timing_ratios.tsv (copy of final_comparison/scores/timing/timing_ratios.tsv)\n')
    fh.write('# note: main-cluster compute node (AMD EPYC 9654), shared with other jobs; ratio = pipeline wall / PLUCK wall in the same paired block. kind=gbrs: one repetition, timing only. Timed PLUCK used the scaled per-sample switch (--switch sample-scaled).\n')
    keep.to_csv(fh, sep='\t', index=False)

# ---------------------------------------------------------------- Table S11
e6 = read('fig6a_e6_extra_animals.tsv')
tl = read('figS10_e6_tail.tsv')
rf = read('figS10_e6_read_factor.tsv')
P('### Table S11. Power model: extra animals, extra reads and the tail, every cohort-tissue\n')
P('| Cohort | Test, n | cis-eQTL or strain-effect detections, PLUCK / Salmon-pan | Extra animals, detection-matched | Noncentrality ratio (strong genes) | Median residual-variance ratio | Genes >2× / >4× residual variance | Salmon-ref genes >2× | δ (median) | γ (median) | Read factor for equal noise |')
P('|---|---|---|---|---|---|---|---|---|---|---|')
for coh, lab in [('Keele kidney', 'Keele kidney, 50 bp SE'), ('Keele liver', 'Keele liver, 50 bp SE†'),
                 ('Munger founders', 'Munger founders, 100 bp SE'), ('GSE199702 liver', 'GSE199702 liver, 101 bp SE'),
                 ('GSE199702 heart', 'GSE199702 heart, 101 bp SE'), ('GSE199702 kidney', 'GSE199702 kidney, 101 bp SE')]:
    a = e6[e6.cohort == coh].iloc[0]
    t = tl[tl.cohort == coh].iloc[0]
    r = rf[rf.source == coh].iloc[0]
    xn = 'n/a (no subsampling)' if pd.isna(a.extra_n_detection_matched) else f'{a.extra_n_detection_matched:.3f}'
    P(f'| {lab} | {a.test}, {a.n_samples} | {a.detections_hybrid:,} / {a.detections_salmon_pan:,} | {xn} | {a.noncentrality_ratio_strong_genes:.3f} | '
      f'{a.median_residual_variance_ratio:.3f} | {100*t.pluck_frac_gt2:.1f}% / {100*t.pluck_frac_gt4:.1f}% | {100*t.salmon_ref_frac_gt2:.1f}% | '
      f'{r.delta_median:.2f} | {r.gamma_median:.5f} | {r.read_factor_text.replace("-", "–").replace("<=", "≤")} |')
P('')
P('Read simulation (independent expected-count truth): δ 1.11 / 0.80 / 0.57 and read factor 1.51 / 1.38 / 1.27 at 50 bp SE / 101 bp SE / 101 bp PE, with no floor (γ ≈ 0).\n')
g = read('figS10_e6_extra_n_grid.tsv')
P('**Animals PLUCK needs relative to Salmon-pan for one gene with T reads per sample (calculator; power model)**\n')
P('| Layout and fit | Residual SD (log2) | T = 50 | 100 | 300 | 1,000 |')
P('|---|---|---|---|---|---|')
for lay, col, name in [('50 SE', 'extra_n_real_median_gene', '50 bp SE, Keele liver real fit'), ('50 SE', 'extra_n_sim', '50 bp SE, simulation'),
                       ('101 SE', 'extra_n_real_median_gene', '101 bp SE, Munger real fit'), ('101 SE', 'extra_n_sim', '101 bp SE, simulation'),
                       ('101 PE', 'extra_n_sim', '101 bp PE, simulation')]:
    for sd in (0.15, 0.25, 0.35):
        s = g[(g.layout == lay) & (g.resid_sd == sd)].sort_values('T')
        P(f'| {name} | {sd:.2f} | ' + ' | '.join(f'{v:.3f}' for v in s[col]) + ' |')
P('')

# ---------------------------------------------------------------- Table S12
sm = read('fig6a_r7_scan_summary.tsv')
avt = read('fig6b_animals_vs_time.tsv')
ex13 = read('fig6a_r7_sensitivity_excl13.tsv')
P('### Table S12. Strains reduced-read Salmon needs to match PLUCK’s cis-eQTL yield (Keele, 50 bp, 47 strains)\n')
P('| Tissue | Arm | Detections at 47 strains | Full-read Salmon-pan strains giving the same detections | Extra animals vs full-read Salmon-pan (MC SE) | Strains relative to PLUCK, proportional reading | Strains needed to match PLUCK, extrapolated (linear / log-log; inference) | Wall time / PLUCK, 4 threads, 50 bp samples |')
P('|---|---|---|---|---|---|---|---|')
ARM2LAB = {'hybrid (frozen v2 core)': 'PLUCK', 'Salmon-pan, 25% of reads': 'Salmon-pan 25%', 'Salmon-pan, 10% of reads': 'Salmon-pan 10%',
           'Salmon-ref, 25% of reads': 'Salmon-ref 25%', 'Salmon-ref, 10% of reads': 'Salmon-ref 10%', 'Salmon-ref, 100%': 'Salmon-ref'}
for tis, tname in [('keele_kidney', 'Keele kidney'), ('keele_liver', 'Keele liver†')]:
    for r in sm[(sm.tissue == tis) & (sm.variant == 'asis')].itertuples():
        lab = ARM2LAB[r.arm]
        if lab == 'Salmon-ref':
            P(f'| {tname} | Salmon-ref (all reads) | {r.detections_at_N:,} | not defined (more detections than Salmon-pan) | not defined | not defined | not defined | 16.2–17.8 |')
            continue
        v = avt[(avt.tissue == tis) & (avt.tool == ('hybrid' if lab == 'PLUCK' else lab))]
        tt = '1' if lab == 'PLUCK' else f'{v.wall_ratio_to_hybrid_4thr_median_min_over_samples.iloc[0]:.1f}–{v.wall_ratio_to_hybrid_4thr_median_max_over_samples.iloc[0]:.1f}'
        ext = '' if lab == 'PLUCK' else (f'{r.strains_needed_to_match_hybrid_extrap_lin:g} / {r.strains_needed_to_match_hybrid_extrap_loglog:g}'
                                          if not pd.isna(r.strains_needed_to_match_hybrid_extrap_lin) else f'not quoted / {r.strains_needed_to_match_hybrid_extrap_loglog:g} (too far beyond 47)')
        P(f'| {tname} | {lab} | {r.detections_at_N:,} | {r.pan_equiv_n:.2f} | {r.extra_n_vs_full_salmon_pan:.3f} ({r.mc_se:.3f}) | '
          f'{r.ratio_to_hybrid_proportional:.3f} | {ext} | {tt} |')
P('')
P('Salmon-pan on all reads: 3,618 (kidney) and 2,360 (liver) detections at 47 strains; wall time 73.8–81.0× PLUCK at 4 threads on the 50 bp samples.\n')
P('**Sensitivity analyses (ordering of the arms unchanged in both)**\n')
P('| Tissue | Arm | Rescaled pseudocount: extra animals vs full-read Salmon-pan | Rescaled: strains relative to PLUCK | Without the 13 low-mapping libraries: N | Extra animals vs full-read Salmon-pan | Strains relative to PLUCK |')
P('|---|---|---|---|---|---|---|')
EXCODE = {'PLUCK': 'hybrid', 'Salmon-pan 25%': 'pan_p25', 'Salmon-pan 10%': 'pan_p10', 'Salmon-ref 25%': 'ref_p25', 'Salmon-ref 10%': 'ref_p10'}
for tis, tname in [('keele_kidney', 'Keele kidney'), ('keele_liver', 'Keele liver†')]:
    for r in sm[(sm.tissue == tis) & (sm.variant == 'rescaled')].itertuples():
        lab = ARM2LAB[r.arm]
        if lab == 'Salmon-ref':
            continue
        e = ex13[(ex13.tissue == tis) & (ex13.arm == EXCODE[lab])].iloc[0]
        P(f'| {tname} | {lab} | {r.extra_n_vs_full_salmon_pan:.3f} | {r.ratio_to_hybrid_proportional:.3f} | {e.N} | {e.extra_n_vs_full_salmon_pan:.3f} | {e.ratio_to_hybrid_proportional:.3f} |')
P('')
print('\n'.join(out))
