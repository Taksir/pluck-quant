#!/usr/bin/env python3
"""Supplementary Fig. S3, the planted cis-eQTL simulation: every check (G1-G5), the null-relative check and the added
analyses. Output figures/FigS3.pdf/.svg. Reads only figures/data/; data file per panel in the comment above it."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fig2_style import P, FOUNDER, FOUNDERS, GREY, MM, DATA, FS_ANN, setup, mstyle, panel_letter, save

setup()
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
LL = {'se50': '50 bp SE', 'se101': '101 bp SE'}
M4 = ['hybrid', 'salmon_pan', 'salmon_ref', 'ref']
M6 = ['hybrid', 'vf', 'emg', 'salmon_pan', 'salmon_ref', 'ref']
fig = plt.figure(figsize=(180 * MM, 218 * MM))
gs = fig.add_gridspec(4, 1, hspace=0.75, left=0.08, right=0.99, top=0.95, bottom=0.05)
hl = [Line2D([], [], ls='none', label=P[c][0], marker=P[c][2], ms=5.2 if c == 'hybrid' else 4.2, mfc=P[c][1] if P[c][3] else 'white',
             mec='white' if c == 'hybrid' else P[c][1], mew=0.5 if c == 'hybrid' else 0.8) for c in M6]
fig.legend(handles=hl, loc='upper center', ncol=6, frameon=False, bbox_to_anchor=(0.53, 0.998), columnspacing=1.4, handletextpad=0.3)
fdr = rd('figS3_idea1_founder.tsv')

# row 1: a power by |beta| (G2; data: fig2c_idea1_power.tsv), b effect-size shift by |beta| (G3; data: figS3_idea1_effect.tsv)
g1 = gs[0].subgridspec(1, 3, width_ratios=[1, 1, 1.3], wspace=0.35)
pw = rd('fig2c_idea1_power.tsv')
B = ['0.25', '0.5', '0.75', '1']
for k, pop in enumerate(['cc', 'do']):
    ax = fig.add_subplot(g1[k])
    for lay, ls in (('se50', '-'), ('se101', '--')):
        r = pw[(pw.population == pop) & (pw.layout == lay)].iloc[0]
        for who, c in (('pan', 'salmon_pan'), ('hybrid', 'hybrid')):
            y = [r[f'power_b{b}_{who}'] for b in B]
            ax.plot([float(b) for b in B], y, color=P[c][1], ls=ls, lw=1.2 if c == 'hybrid' else 0.8)
            ax.scatter([float(b) for b in B], y, **mstyle(c, 8))
        tx, ty, ha = (0.97, 0.13, 'right') if pop == 'cc' else (0.03, 0.84, 'left')
        ax.text(tx, ty - (0.09 if lay == 'se101' else 0), f'{LL[lay].replace(" SE", "")}: {r.hybrid_detected} vs {r.pan_detected} '
                f'({r.ratio:.3f})', transform=ax.transAxes, fontsize=FS_ANN, ha=ha)
    hdr = ['detections,', 'PLUCK vs Salmon-pan'] if pop == 'cc' else ['detections, PLUCK vs Salmon-pan']   # CC: two lines clear the curve
    for i, h in enumerate(reversed(hdr)):
        ax.text(tx, ty + 0.09 * (i + 1), h, transform=ax.transAxes, fontsize=FS_ANN, ha=ha)
    ax.set_xlabel('Planted |β| (log$_2$)'); ax.set_ylim(0, 1.05); ax.set_xticks([0.25, 0.5, 0.75, 1])
    ax.set_title(('Detections: CC, 64 animals' if pop == 'cc' else 'Detections: DO, 32 animals'), fontsize=7)
    if k == 0:
        ax.set_ylabel('Power at BH 5%'); panel_letter(ax, 'a', x=-0.3)
    else:   # line-style key in the DO facet, where the curves stay low
        ax.legend(handles=[Line2D([], [], color='#444444', ls='-', lw=0.8, label='50 bp SE'),
                           Line2D([], [], color='#444444', ls='--', lw=0.8, label='101 bp SE')], loc='center left', bbox_to_anchor=(0.0, 0.5), frameon=False, fontsize=FS_ANN)
ax = fig.add_subplot(g1[2])
ef = rd('figS3_idea1_effect.tsv'); ef = ef[(ef.source == 'locked') & ef.beta.notna()]
off = dict(zip(M6, np.linspace(-0.05, 0.05, 6)))
ax.axhline(0, color='black', lw=0.5)
for lay, ls in (('se50', '-'), ('se101', '--')):
    for c in M6[::-1]:
        d = ef[(ef.layout == lay) & (ef.method == c)].sort_values('beta')
        x = d.beta + off[c] + (0.012 if lay == 'se101' else -0.012)
        kw = mstyle(c, 8)
        if lay == 'se101':
            kw.update(alpha=0.6)
        ax.scatter(x, d.median_shift, **kw)
ax.set_ylim(-0.06, 0.06); ax.set_xticks([0.25, 0.5, 0.75, 1]); ax.set_xlabel('Planted |β| (log$_2$)')
ax.set_ylabel('Median effect-size shift (log$_2$)'); ax.set_title('Effect-size shift, CC (faded: 101 bp)', fontsize=7)
panel_letter(ax, 'b', x=-0.22)

# row 2: c per-founder median shift (G4a) and d share of genes beyond 0.05 log2 (G4b), CC (data: figS3_idea1_founder.tsv, original-calibration rows)
g2 = gs[1].subgridspec(1, 5, width_ratios=[1, 1, 0.28, 1, 1], wspace=0.12)
loc = fdr[fdr.source == 'locked']
o4 = dict(zip(M4, [-0.27, -0.09, 0.09, 0.27]))
axc, axd = [], []
for j, lay in enumerate(['se50', 'se101']):
    ax = fig.add_subplot(g2[j], sharey=axc[0] if axc else None); axc.append(ax)
    ax.axhline(0, color='black', lw=0.5)
    for c in M4[::-1]:
        d = loc[(loc.layout == lay) & (loc.method == c)].set_index('founder').loc[FOUNDERS]
        ax.scatter(np.arange(7) + o4[c], d.median_delta, **mstyle(c, 10))
    ax.set_xticks(range(7)); ax.set_xticklabels([FOUNDER[f] for f in FOUNDERS], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_title(f'Founder shift,\nCC, {LL[lay]}', fontsize=7)
    if j:
        plt.setp(ax.get_yticklabels(), visible=False)
    else:
        ax.set_ylabel('Median founder shift (log$_2$)'); panel_letter(ax, 'c', x=-0.4)
axc[0].set_ylim(-0.2, 0.03)
for j, lay in enumerate(['se50', 'se101']):
    ax = fig.add_subplot(g2[j + 3], sharey=axd[0] if axd else None); axd.append(ax)
    for c in M4[::-1]:
        d = loc[(loc.layout == lay) & (loc.method == c)].set_index('founder').loc[FOUNDERS]
        ax.scatter(np.arange(7) + o4[c], 100 * d.frac_beyond_bound, **mstyle(c, 10))
    ax.set_yscale('log'); ax.set_ylim(0.5, 80); ax.set_yticks([1, 10, 50]); ax.set_yticklabels(['1%', '10%', '50%'])
    ax.set_xticks(range(7)); ax.set_xticklabels([FOUNDER[f] for f in FOUNDERS], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_title(f'Hit genes beyond 0.05 log$_2$\nCC, {LL[lay]}', fontsize=7)
    if j:
        plt.setp(ax.get_yticklabels(), visible=False)
    else:
        ax.set_ylabel('Hit genes beyond 0.05 log$_2$'); panel_letter(ax, 'd', x=-0.4)

# row 3: e null-relative check (data: figS3_idea1_founder.tsv, reported only), f G4c DO pairs (data: figS3_idea1_pairs.tsv)
g3 = gs[2].subgridspec(1, 3, width_ratios=[1, 1, 0.8], wspace=0.38)
nr = fdr[fdr.source.str.startswith('null-relative')]
axe = []
for j, lay in enumerate(['se50', 'se101']):
    ax = fig.add_subplot(g3[j], sharey=axe[0] if axe else None); axe.append(ax)
    ax.axhline(0, color='black', lw=0.5)
    for c in ['salmon_ref', 'salmon_pan', 'hybrid']:
        d = nr[(nr.layout == lay) & (nr.method == c)].set_index('founder').loc[FOUNDERS]
        x = np.arange(7) + o4[c]
        ax.errorbar(x, 100 * d.excess, yerr=200 * d.excess_se, fmt='none', ecolor=P[c][1], elinewidth=0.6)
        ax.scatter(x, 100 * d.excess, **mstyle(c, 10))
    ax.set_xticks(range(7)); ax.set_xticklabels([FOUNDER[f] for f in FOUNDERS], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_title(f'CC, {LL[lay]}: excess over own null (±2 SE)', fontsize=7)
    if j:
        plt.setp(ax.get_yticklabels(), visible=False)
    else:
        ax.set_ylabel('Share beyond bound minus\nB6-split null (percentage points)'); panel_letter(ax, 'e', x=-0.3)
axe[0].set_ylim(-3, 40)
ax = fig.add_subplot(g3[2])
pr = rd('figS3_idea1_pairs.tsv'); pr = pr[pr.source == 'locked']
rng = np.random.default_rng(4)
ax.axhline(0, color='black', lw=0.5)
for j, lay in enumerate(['se50', 'se101']):
    v = pr[pr.layout == lay].median_delta.values
    ax.scatter(j + (rng.random(len(v)) - 0.5) * 0.4, v, **mstyle('hybrid', 8))
    ax.text(j, 0.0165, f'max |median|\n{np.abs(v).max():.4f}', fontsize=FS_ANN, ha='center')
ax.set_xticks([0, 1]); ax.set_xticklabels(['50 bp SE', '101 bp SE']); ax.set_xlim(-0.6, 1.6); ax.set_ylim(-0.025, 0.025)
ax.set_ylabel('PLUCK median shift (log$_2$)'); ax.set_title('Heterozygous founder\npairs, DO (28)', fontsize=7)
panel_letter(ax, 'f', x=-0.42)

# row 4: g G5 call rate by branch (data: figS3_idea1_branch.tsv), h context statistic, original calibration vs recalibration
# (data: figS3_idea1_effect.tsv M2_null_cc rows) and i the fixed-unit call-rate plane at 101 bp (data: fig2d_g5_fixed_unit.tsv)
g4 = gs[3].subgridspec(1, 3, width_ratios=[1.25, 0.9, 1.0], wspace=0.42)
ax = fig.add_subplot(g4[0])
br = rd('figS3_idea1_branch.tsv'); br = br[br.branch.isin(['gene_probe_branch', 'EM_branch'])]
cats = [(lay, sw, b) for lay in ('se50', 'se101') for sw in ('swM', 'swR') for b in ('gene_probe_branch', 'EM_branch')]
for i, (lay, sw, b) in enumerate(cats):
    for c, dx in (('hybrid', -0.15), ('salmon_pan', 0.15)):
        r = br[(br.layout == lay) & (br.switch == sw) & (br.branch == b) & (br.method == c)].iloc[0]
        ax.scatter([i + dx], [100 * r.call_rate], **mstyle(c, 10))
ax.set_xticks(range(len(cats)))
ax.set_xticklabels([f"{'copy-pres.' if sw == 'swM' else 'read-pres.'}, {'probe' if b.startswith('gene') else 'EM'}" for _, sw, b in cats],
                   rotation=45, ha='right', rotation_mode='anchor')
ax.axvline(3.5, color=GREY, lw=0.4)
ax.text(1.5, 58, '50 bp SE', ha='center', fontsize=FS_ANN); ax.text(5.5, 58, '101 bp SE', ha='center', fontsize=FS_ANN)
ax.set_ylim(-2, 62); ax.set_ylabel('Switch genes called (%)')
ax.set_title('Isoform switches by branch', fontsize=7)
panel_letter(ax, 'g', x=-0.25)

ax = fig.add_subplot(g4[1])
e2 = rd('fig2a_e2_metrics_all.tsv')
ctx = rd('figS3_idea1_effect.tsv'); ctx = ctx[ctx.M2_null_cc.notna()]
for j, lay in enumerate(['se50', 'se101']):
    for c, dx in (('hybrid', -0.15), ('salmon_pan', 0.15)):
        a = ctx[(ctx.layout == lay) & (ctx.method == c) & (ctx.source == 'locked')].M2_null_cc.iloc[0]
        b = ctx[(ctx.layout == lay) & (ctx.method == c) & (ctx.source != 'locked')].M2_null_cc.iloc[0]
        ax.annotate('', xy=(j + dx, b), xytext=(j + dx, a), arrowprops=dict(arrowstyle='-|>', color=P[c][1], lw=0.6, mutation_scale=5))
        ax.scatter([j + dx], [a], **mstyle(c, 10))
        ref = e2[(e2.cond == 'clean') & (e2.layout == lay) & (e2.truth == 'expected') & (e2.tier == 'core') & (e2.method == c)].M2_err_B6.iloc[0]
        ax.plot([j + dx - 0.1, j + dx + 0.1], [ref, ref], color='black', lw=0.9)
ax.set_xticks([0, 1]); ax.set_xticklabels(['50 bp SE', '101 bp SE']); ax.set_xlim(-0.5, 1.5); ax.set_ylim(0.03, 0.12)
ax.set_ylabel('Median |log$_2$ error|, null genes')
ax.set_title('Original → recalibrated\n(black bar: simulation value)', fontsize=7)
panel_letter(ax, 'h', x=-0.42)

ax = fig.add_subplot(g4[2])
g5 = rd('fig2d_g5_fixed_unit.tsv'); cr = g5[(g5.layout == 'se101') & g5.analysis.str.startswith('call rate')]
for unit, code, sty, lab in [('tpm', 'salmon_ref', 'salmon_ref', 'Salmon-ref, TPM'), ('tpm', 'salmon_pan', 'salmon_pan', 'Salmon-pan, TPM'),
                             ('reads', 'vf', 'vf', 'gene probes alone'), ('reads', 'ref', 'ref', 'ref-unique probes'),
                             ('reads', 'emg', 'emg', 'EM alone'), ('reads', 'salmon_ref', 'salmon_ref', 'Salmon-ref, reads'),
                             ('reads', 'salmon_pan', 'salmon_pan', 'Salmon-pan, reads'), ('reads', 'hybrid', 'hybrid', 'PLUCK')]:
    r = cr[(cr.unit == unit) & (cr.tool_code == code)].iloc[0]
    ax.scatter([100 * r.swM_call_rate], [100 * r.swR_call_rate], **mstyle(sty, 10))
t = cr[(cr.unit == 'reads') & (cr.tool_code == 'truth')].iloc[0]
ax.scatter([100 * t.swM_call_rate], [100 * t.swR_call_rate], marker='*', s=20, c='black', zorder=6)
ax.text(27.5, 12, 'PLUCK', color=P['hybrid'][1], fontsize=FS_ANN, ha='right', va='center')
ax.text(4, 25, 'Salmon in TPM', fontsize=FS_ANN, ha='left', va='top', color=GREY)
ax.text(39.5, 15, 'Salmon in\nread counts', fontsize=FS_ANN, ha='right', va='bottom', color=GREY)
ax.set_xlim(-1.5, 40); ax.set_ylim(-1.5, 40)
ax.set_xlabel('Calls, copy-pres. switches (%)'); ax.set_ylabel('Calls, read-pres. switches (%)')
ax.set_title('Fixed units, CC 101 bp', fontsize=7)
panel_letter(ax, 'i', x=-0.3)
save(fig, 'FigS3')
print('saved FigS3')
