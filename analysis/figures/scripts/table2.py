#!/usr/bin/env python3
"""Table 2: headline grid, pipeline x cohort accuracy (part A) and per-sample cost at four threads (part B).
Reads only figures/data/:
  table2_accuracy_frameA.tsv, table2_accuracy_rnaskim_frame.tsv     (part A)
  table2_timing_4thr.tsv, fig6c_rnaskim_timing_second_cluster.tsv, fig6c_rnaskim_peak_rss.tsv, fig6c_timing_ratios.tsv (part B)
Writes tables/Table2.tsv (both parts, machine-readable) and prints the same table as Markdown.
Run: python3 -I figures/scripts/table2.py"""
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
D = '/mnt/project-files/paper/figures/data'
OUT = '/mnt/project-files/paper/tables/Table2.tsv'
ORDER = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags', 'RNA-Skim',
         'ref-unique probes', 'GBRS']
COH = ['Keele kidney', 'Keele liver', 'Munger founders']
CL = {'Keele kidney': 'Keele kidney, 50 bp', 'Keele liver': 'Keele liver, 50 bp†', 'Munger founders': 'Munger founders, 100 bp'}


def rd(n):
    return pd.read_csv(f'{D}/{n}', sep='\t', comment='#')


lab = lambda t: 'PLUCK' if t == 'hybrid' else t
a = rd('table2_accuracy_frameA.tsv'); a['tool'] = a.tool.map(lab); a['frame_note'] = ''
r = rd('table2_accuracy_rnaskim_frame.tsv').rename(columns={'cohort_label': 'cohort'})
r = r[r.tool == 'RNA-Skim'].copy(); r['frame_note'] = 'own gene frame (P_rnaskim)'
acc = pd.concat([a, r])


def ci(s):
    lo, hi = [float(x) for x in s.strip('[]').split(',')]
    return lo, hi


def sgn(x, d=4):
    return ('+' if x > 0 else ('−' if x < 0 else '')) + f'{abs(x):.{d}f}'


rows, md = [], []
md.append('| Pipeline | Within-gene Spearman | Noise SD | Median Δ [95% CI] | Spurious rate (%) |')
md.append('|---|---|---|---|---|')
for c in COH:
    ng = int(acc[(acc.cohort == c) & (acc.tool == 'PLUCK')].n_delta_genes.iloc[0])
    md.append(f'| *{CL[c]}: 47 samples, {ng:,} genes in Δ and spurious rate* | | | | |')
    for t in ORDER:
        x = acc[(acc.cohort == c) & (acc.tool == t)]
        if t == 'GBRS':
            continue
        if x.empty:
            md.append(f'| {t} | not run | | | |')
            rows.append(dict(part='A', pipeline=t, cohort=CL[c], status='not run'))
            continue
        x = x.iloc[0]; lo, hi = ci(x.delta_ci95)
        mark = '§' if t == 'RNA-Skim' else ''
        name = f'**{t}**' if t == 'PLUCK' else t
        md.append(f'| {name}{mark} | {x.within_gene_spearman:.4f} | {x.noise_sd:.4f} | {sgn(x.median_delta)} [{sgn(lo)}, {sgn(hi)}] | {100 * x.spurious_rate:.1f} |')
        rows.append(dict(part='A', pipeline=t, cohort=CL[c], within_gene_spearman=x.within_gene_spearman, noise_sd=x.noise_sd,
                         median_delta=x.median_delta, delta_ci95=x.delta_ci95, spurious_rate=x.spurious_rate, n_delta_genes=x.n_delta_genes,
                         frame='P_rnaskim' if t == 'RNA-Skim' else 'A'))

# part B: cost per sample at 4 threads
t4 = rd('table2_timing_4thr.tsv'); t4['tool'] = t4.tool.map(lab)
sy = rd('fig6c_rnaskim_timing_second_cluster.tsv'); rss = rd('fig6c_rnaskim_peak_rss.tsv')
gb = rd('fig6c_timing_ratios.tsv'); gb = gb[(gb.kind == 'gbrs') & (gb.method == 'gbrs') & (gb['T'] == 4)]
SAMPLES = [('50', 'CC021'), ('50', 'CC001'), ('101', 'LiverCC001F')]
SL = {'CC021': 'Keele liver CC021, 50 bp', 'CC001': 'Keele liver CC001, 50 bp', 'LiverCC001F': 'GSE199702 liver CC001 F, 101 bp'}
md2 = ['| Pipeline | ' + ' | '.join(f'{SL[s]}: wall (s), ratio' for _, s in SAMPLES) + ' | Peak RAM (GB) |',
       '|---|' + '---|' * (len(SAMPLES) + 1)]
for t in ORDER:
    cells, rams = [], []
    for lay, s in SAMPLES:
        if t == 'RNA-Skim':
            x = sy[(sy['sample'] == s) & (sy.other == 'hybrid')]
            if x.empty:
                cells.append('not timed'); continue
            x = x.iloc[0]
            cells.append(f'{x.median_wall_rnaskim_s:.0f} (PLUCK {x.median_wall_other_s:.1f}), {x.median_ratio_rnaskim_over_other:.0f}×')
            rams.append(float(rss[(rss['sample'] == s) & (rss.method == 'rnaskim')].rss_mb.iloc[0]))
            rows.append(dict(part='B', pipeline=t, sample=SL[s], machine='second compute cluster', wall_median_s=x.median_wall_rnaskim_s, pluck_wall_median_s_same_block=x.median_wall_other_s,
                             ratio_median=x.median_ratio_rnaskim_over_other, ratio_min=x.min_ratio, ratio_max=x.max_ratio, n_reps=x.n_reps))
            continue
        if t == 'GBRS':
            x = gb[gb['sample'] == s]
            if x.empty:
                cells.append('not timed'); continue
            x = x.iloc[0]
            cells.append(f'{x.wall_median_s:.1f}, {x.ratio_median:.1f}× (1 rep)')
            rams.append(float(x.rss_mb_median))
            rows.append(dict(part='B', pipeline=t, sample=SL[s], machine='compute node', wall_median_s=x.wall_median_s, ratio_median=x.ratio_median,
                             n_reps=1, rss_mb_median=x.rss_mb_median))
            continue
        x = t4[(t4['sample'] == s) & (t4.tool == t)]
        if x.empty:
            cells.append('not timed' if t != 'kallisto-pan' else 'skipped'); continue
        x = x.iloc[0]
        cells.append(f'{x.wall_median_s:.1f}, 1' if t == 'PLUCK' else f'{x.wall_median_s:.1f}, {x.ratio_median:.1f}×')
        rams.append(float(x.rss_mb_median))
        rows.append(dict(part='B', pipeline=t, sample=SL[s], machine='compute node', wall_median_s=x.wall_median_s, ratio_median=x.ratio_median,
                         ratio_min=x.ratio_min, ratio_max=x.ratio_max, cpu_ratio_median=x.cpu_ratio_median, rss_mb_median=x.rss_mb_median, n_reps=x.n_reps))
    if t == 'ref-unique probes':
        cells = ['not timed (emulation)'] * len(SAMPLES)
    ram = '' if not rams else (f'{min(rams) / 1000:.1f}' if round(min(rams) / 1000, 1) == round(max(rams) / 1000, 1)
                                else f'{min(rams) / 1000:.1f}–{max(rams) / 1000:.1f}')
    name = f'**{t}**' if t == 'PLUCK' else (t + '‡' if t == 'RNA-Skim' else t)
    md2.append(f'| {name} | ' + ' | '.join(cells) + f' | {ram} |')

pd.DataFrame(rows).to_csv(OUT, sep='\t', index=False)
print('\n'.join(md)); print(); print('\n'.join(md2))
