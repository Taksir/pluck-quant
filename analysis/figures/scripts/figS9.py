#!/usr/bin/env python3
"""Fig. S9: structural-variant check of the probe panel against the seven non-B6 founder assemblies.
Reads only figures/data/. Writes figures/FigS9.pdf and FigS9.svg."""
import numpy as np, pandas as pd
from u5_style import *

setup()
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
fig = plt.figure(figsize=(180 * MM, 127 * MM))
ORANGE = PIPE['PLUCK'][0]
FOUNDERS = [('A', 'A/J'), ('C', '129S1'), ('D', 'NOD'), ('E', 'NZO'), ('F', 'CAST'), ('G', 'PWK'), ('H', 'WSB')]

# a: data figS9a_sv_probe_status_u5.tsv
seg = rd('figS9a_sv_probe_status_u5.tsv').set_index('segment')
SEGS = [('core50_gene', 'core gene probes, 50 bp', 'black', 'o', True), ('core101_gene', 'core gene probes, 101 bp', 'black', 'o', False),
        ('core50_EM', 'core EM probes, 50 bp', '#7A7A7A', 's', True), ('ext101', 'extended tier, 101 bp', '#7A7A7A', 's', False),
        ('ctrl_ref_only_variant_probes', 'control: probes carrying an annotated variant', '#555555', 'x', True)]
ax = ax_mm(fig, 16, 9, 74, 48)
xs = np.arange(len(FOUNDERS) + 1)
for k, (sg, lab, col, mk, filled) in enumerate(SEGS):
    vals = [seg.loc[sg, f'{f}_absent'] * 100 for f, _ in FOUNDERS] + [seg.loc[sg, 'N_absent'] * 100]
    dx = (k - 2) * 0.1
    ax.plot(xs + dx, vals, ls='none', marker=mk, ms=3.6, mfc=col if filled else 'white', mec=col, mew=0.8, label=lab)
    floor = seg.loc[sg, 'N_absent'] * 100
ax.set_yscale('log'); ax.set_ylim(0.25, 25); ax.set_yticks([0.3, 1, 3, 10]); ax.set_yticklabels(['0.3%', '1%', '3%', '10%']); ax.minorticks_off()
ax.set_xticks(xs); ax.set_xticklabels([n for _, n in FOUNDERS] + ['B6NJ\n(floor)'], linespacing=1.0)
ax.axvline(len(FOUNDERS) - 0.5, color='#CCCCCC', lw=0.6)
ax.set_ylabel('Single-copy probes absent\nfrom the founder assembly', linespacing=1.1)
ax.legend(loc='upper left', frameon=False, handletextpad=0.3, borderaxespad=0.2, labelspacing=0.25, fontsize=7, ncol=1)
letter(fig, 1, 5, 'a')
mtext(fig, 180 - 3, 3, 'No other analysis uses this measurement', ha='right', va='top', color='#333333', style='italic')

# b: data figS9b_sv_vs_e5_u5.tsv (threshold 0.25)
sve = rd('figS9b_sv_vs_e5_u5.tsv'); sve = sve[sve.th == 0.25]
ax = ax_mm(fig, 110, 9, 66, 48)
for i, (f, n) in enumerate(FOUNDERS):
    for m, lab, dx in (('vf', 'PLUCK', -0.14), ('salmon_pan', 'Salmon-pan', 0.14)):
        r = sve[(sve.founder == f) & (sve.method == m)].iloc[0]
        ax.plot([i + dx, i + dx], [r.median_abs_other, r.median_abs_flagged], color='#AAAAAA', lw=0.7, zorder=2)
        pmark(ax, i + dx, r.median_abs_other, lab, size=3.4, scale=0.75, alpha=0.55)   # all other genes: same marker, smaller and lighter
        pmark(ax, i + dx, r.median_abs_flagged, lab, size=3.4)                          # SV-affected genes: full size
    ax.text(i, 0.415, str(int(sve[(sve.founder == f) & (sve.method == 'vf')].n_flagged_genes.iloc[0])), ha='center', va='bottom', color='#333333')
ax.text(-0.8, 0.415, 'genes:', ha='right', va='bottom', color='#333333')
ax.set_xticks(range(7)); ax.set_xticklabels([n for _, n in FOUNDERS]); ax.set_xlim(-0.6, 6.6)
ax.set_ylim(0, 0.45); ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
ax.set_ylabel('Median |log$_2$ error| on founder mice\n(Munger founders, 100 bp)', linespacing=1.1)
lx = 110
for k, lab in enumerate(['PLUCK', 'Salmon-pan']):
    lax = ax_mm(fig, lx + k * 22, 64.5, 4, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    pmark(lax, 0.5, 0.5, lab, size=3.4)
    mtext(fig, lx + k * 22 + 5, 66, lab, va='center')
mtext(fig, lx, 70, 'full size: SV-affected genes in that founder;\nsmaller, lighter: all other genes', va='center', linespacing=1.1)
letter(fig, 97, 5, 'b')

# c: data tableS14_sv_risk_genes_u5.tsv (number of founders with >= 25% excess loss per gene)
sv = rd('tableS14_sv_risk_genes_u5.tsv')
nf = sv.founders_ge25.fillna('').astype(str).str.len()
cnt = nf.value_counts().sort_index()
ax = ax_mm(fig, 16, 78, 52, 38)
ax.bar(cnt.index, cnt.values, width=0.7, color='#7A7A7A')
for x, v in cnt.items():
    ax.text(x, v + 6, str(v), ha='center', va='bottom')
ax.set_xticks(range(1, 8)); ax.set_xlim(0.4, 7.6); ax.set_ylim(0, cnt.max() * 1.18)
ax.set_xlabel('Founders in which the gene loses ≥25% of probes'); ax.set_ylabel('SV-risk core genes')
letter(fig, 1, 74, 'c')

# d: data figS9c_sv_flag_overlap_u5.tsv
ov = rd('figS9c_sv_flag_overlap_u5.tsv').set_index('set')
ax = ax_mm(fig, 96, 80, 70, 30)
rows = [('SV-risk genes (>= 25%)', 'SV-risk genes (n = 567)'), ('all core genes', 'all core genes (n = 16,885)')]
for j, (key, lab) in enumerate(rows):
    for k, (col, l2, face) in enumerate((('flagged_50bp_frac', '50 bp design', ORANGE), ('flagged_100bp_frac', '100/101 bp design', 'white'))):
        y = j * 1.0 + (-0.18 if k == 0 else 0.18)
        v = ov.loc[key, col] * 100
        ax.barh(y, v, height=0.32, color=face, edgecolor=ORANGE, lw=0.8)
        ax.text(v + 1.5, y, '%.1f%%' % v, va='center')
ax.set_yticks([0, 1]); ax.set_yticklabels([r[1] for r in rows]); ax.tick_params(axis='y', length=0)
ax.set_ylim(1.55, -0.55); ax.set_xlim(0, 100); ax.set_xticks([0, 25, 50, 75, 100]); ax.set_xticklabels(['0', '25', '50', '75', '100%'])
ax.set_xlabel('Share flagged by the reliability flag ($u_g$ < 0.6)')
for k, (l2, face) in enumerate((('50 bp design', ORANGE), ('100/101 bp design', 'white'))):
    lax = ax_mm(fig, 120 + k * 28, 120.5, 4, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    lax.add_patch(plt.Rectangle((0.1, 0.15), 0.8, 0.7, fc=face, ec=ORANGE, lw=0.8))
    mtext(fig, 120 + k * 28 + 5, 122, l2, va='center')
letter(fig, 76, 74, 'd')

save(fig, 'FigS9')
print('wrote FigS9')
