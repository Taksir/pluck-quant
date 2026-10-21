#!/usr/bin/env python3
"""Fig. 1: probe design, counting, the hybrid switch, panel coverage and the measurement-model check.
Reads only figures/data/ (written by fig1_data.py) and writes figures/Fig1.pdf and Fig1.svg.
Colours: PLUCK #D55E00; ref-unique probes #332288; gene-probe branch #F08A5D and EM branch #8C2D04 (the ablation
hues of Fig. 5e); sans-serif, bold lowercase panel letters 10 pt.
Elsevier artwork rule: all lettering >= 7 pt and sub/superscripts >= 6 pt at the printed size, which is
the size drawn here (180 mm wide). Mathtext scripts are drawn at 0.87 of the base size (7 pt -> 6.1 pt).
Check with: svg_fontsizes.py ../Fig1.svg
usage: fig1.py"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib._mathtext as _mathtext
_mathtext.SHRINK_FACTOR = 0.87   # sub/superscripts 6.1 pt on 7 pt text (default 0.7 gives 4.9 pt)
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon
from matplotlib.lines import Line2D

# Folder constants below point to the original working folder; set them for your copy.
D = '/mnt/project-files/paper/figures/data'
OUT = '/mnt/project-files/paper/figures'
PREVIEW = 'Fig1_preview.png'
MM = 1 / 25.4
FS = 7          # every label; panel letters 10 pt; column headings of panel b 7.5 pt bold
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': FS, 'axes.labelsize': FS, 'axes.titlesize': FS, 'xtick.labelsize': FS, 'ytick.labelsize': FS,
                     'legend.fontsize': FS, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
                     'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic',
                     'mathtext.bf': 'Liberation Sans:bold', 'mathtext.sf': 'Liberation Sans', 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'axes.spines.top': False, 'axes.spines.right': False})
PLUCK = '#D55E00'; REFU = '#332288'; GP = '#F08A5D'; EM = '#8C2D04'
INK = '#222222'; MID = '#7F7F7F'; LIGHT = '#D9D9D9'; PALE = '#F2F2F2'
K = 25
FOUNDER_ROWS = [('A', 'A/J (A)'), ('B', 'C57BL/6J (B, reference)'), ('C', '129S1/SvImJ (C)'), ('D', 'NOD/ShiLtJ (D)'),
                ('E', 'NZO/HlLtJ (E)'), ('F', 'CAST/EiJ (F)'), ('G', 'PWK/PhJ (G)'), ('H', 'WSB/EiJ (H)')]
W_MM, H_MM = 180, 210
# panel rows in mm from the bottom of the figure
ROW_A = (132, 77); ROW_B = (72, 56); ROW_CD = (12, 52)


def read(name):
    return pd.read_csv(f'{D}/{name}', sep='\t', comment='#')


def letter(fig, s, x_mm, y_mm):
    fig.text(x_mm / W_MM, y_mm / H_MM, s, fontsize=10, fontweight='bold', va='top', ha='left')


# ------------------------------------------------------------------------------------------------ panel a
def panel_a(ax, fig):
    # data: fig1a_gene_positions.tsv, fig1a_gene_variants.tsv, fig1a_gene_isoforms.tsv
    # y units are mm; x units are nucleotides of the longer isoform, placed so that x = 0 sits 39 mm from the left edge
    P = read('fig1a_gene_positions.tsv'); V = read('fig1a_gene_variants.tsv'); I = read('fig1a_gene_isoforms.tsv')
    L = int(I.loc[~I.transcript.str.startswith('exon:'), 'length'].max())
    R = 101; W = R - K + 1
    x0_mm, gene_mm = 39.0, 111.0
    upm = L / gene_mm                                       # nucleotides per mm
    ax.set_xlim(-x0_mm * upm, (W_MM - x0_mm) * upm); ax.set_ylim(0, ROW_A[1]); ax.axis('off')
    lab = -1.2 * upm                                        # right edge of the row labels
    ax.text(L / 2, ROW_A[1] - 0.5, '$\\it{Slc25a47}$, 101 bp design: probes at least $R-k+1$ = 77 nt apart, '
            'on the coordinates of the longer TSL1 isoform', ha='center', va='top', fontsize=FS)
    # isoform models: reference positions each TSL1 isoform contains
    isos = I[~I.transcript.str.startswith('exon:')].sort_values('length', ascending=False)
    exons = I[I.transcript.str.startswith('exon:')]
    yi = [69.0, 65.4]
    for k, (tname, grp) in enumerate(isos.groupby('transcript', sort=False)):
        yy = yi[k]
        ax.plot([0, L], [yy, yy], color=LIGHT, lw=0.5, zorder=1)
        for _, r in grp.iterrows():
            ax.add_patch(Rectangle((r.start, yy - 1.1), r.end - r.start, 2.2, fc=MID, ec='none', zorder=2))
        ax.text(lab, yy, f'TSL1 isoform, {int(grp.length.iloc[0]):,} nt', ha='right', va='center', fontsize=FS)
    for _, r in exons.iloc[1:].iterrows():
        ax.plot([r.start, r.start], [yi[0] + 1.3, yi[-1] - 1.3], color='white', lw=0.6, zorder=3)
    # founder haplotypes with annotated SNPs and small indels
    y0, dy = 59.6, 3.55
    for i, (lt, nm) in enumerate(FOUNDER_ROWS):
        yy = y0 - i * dy
        ax.plot([0, L], [yy, yy], color=LIGHT, lw=0.8, zorder=1)
        ax.text(lab, yy, nm, ha='right', va='center', fontsize=FS, color=INK)
        for _, r in V[V.founder == lt].iterrows():
            if r.kind == 'SNP':
                ax.plot([r.tx_start + 0.5] * 2, [yy - 1.2, yy + 1.2], color=INK, lw=0.7, zorder=3, solid_capstyle='butt')
            else:
                ax.add_patch(Rectangle((r.tx_start, yy - 1.2), max(r.tx_end - r.tx_start + 1, 4), 2.4, fc='white', ec=INK, lw=0.6,
                                       zorder=3))
    # candidate 25-mer starts and the first rule each one fails
    order = [('end margin', PALE, 'near an end'), ('not in every isoform', LIGHT, 'not in every isoform'),
             ('shared with another gene', '#A6A6A6', 'in another gene'), ('more than one genome copy', '#595959', '>1 genome copy'),
             ('overlaps a founder variant', REFU, 'on a founder SNP/indel'), ('eligible', 'white', 'eligible')]
    col = {k: c for k, c, _ in order}
    st = P.status.to_numpy(); runs = np.flatnonzero(np.r_[True, st[1:] != st[:-1]]); ends = np.r_[runs[1:], len(st)]
    yc = 27.6
    for s0, e0 in zip(runs, ends):
        ax.add_patch(Rectangle((s0, yc - 1.2), e0 - s0, 2.4, fc=col[st[s0]], ec='none', zorder=2,
                               alpha=0.35 if st[s0] == 'overlaps a founder variant' else 1))
    ax.add_patch(Rectangle((0, yc - 1.2), len(st), 2.4, fc='none', ec=MID, lw=0.4, zorder=3))
    ax.text(lab, yc, '25-mer starts', ha='right', va='center', fontsize=FS)
    # key: only the rules some 25-mer start of this gene fails (none here is in another gene or has >1 genome copy)
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    lx, ly = 0.0, yc - 4.0
    for k, c, txt in order:
        if k not in set(st):
            continue
        ax.add_patch(Rectangle((lx, ly - 1.1), 2.6 * upm, 2.2, fc=c, ec=MID, lw=0.3, alpha=0.35 if k.startswith('overlaps') else 1))
        t = ax.text(lx + 3.6 * upm, ly, txt, va='center', fontsize=FS)
        bb = t.get_window_extent(renderer=rend)
        lx = inv.transform((bb.x1, bb.y0))[0] + 4.5 * upm
    # ref-unique probes (no variant filter), greedy spacing >= R - k + 1; letters = founders whose variants each overlaps
    yr = 14.4
    ref = P[P.ref_unique_probe == 1]
    last, lift = -1e9, 0.0
    for _, r in ref.iterrows():
        hit = r.founder_mask > 0
        ax.add_patch(Rectangle((r.pos, yr - 1.1), K, 2.2, fc='white' if hit else REFU, ec=REFU, lw=0.7, zorder=3))
        if hit:
            lift = 2.6 if (r.pos - last < 110 and lift == 0) else 0.0   # stagger the letters of neighbouring hit probes
            ax.text(r.pos + K / 2, yr + 1.5 + lift, ''.join(sorted(str(r.founders))), ha='center', va='bottom', fontsize=FS,
                    color=REFU)
            last = r.pos
    ax.text(lab, yr, 'ref-unique probes', ha='right', va='center', fontsize=FS, color=REFU)
    nhit = int((ref.founder_mask > 0).sum())
    ax.text(L + 1.5 * upm, yr + 1.0, f'{len(ref)} probes, {nhit} on a\nfounder variant (open)', ha='left', va='center', fontsize=FS,
            color=REFU, linespacing=1.15)
    # PLUCK probes, with the read-start window each one owns
    yp = 7.6
    pl = P[P.pluck_probe == 1]
    for k, (_, r) in enumerate(pl.iterrows()):
        ax.add_patch(Rectangle((r.pos - (R - K), yp - 3.4), W, 1.1, fc=PLUCK if k % 2 == 0 else GP, ec='none', alpha=0.35, zorder=2))
        ax.add_patch(Rectangle((r.pos, yp - 1.1), K, 2.2, fc=PLUCK, ec='white', lw=0.4, zorder=4))
    ax.text(lab, yp, 'PLUCK probes', ha='right', va='center', fontsize=FS, color=PLUCK, fontweight='bold')
    ax.text(lab, yp - 3.0, f'read starts, {W} each', ha='right', va='center', fontsize=FS, color=PLUCK)
    cap = len(pl) * W / np.mean([max(l - R + 1, K) for l in isos.groupby('transcript').length.first()])
    lines(ax, L + 1.5 * upm, yp + 0.6, [f'{len(pl)} probes, none on an', 'annotated SNP/indel;', f'$\\kappa_g$ = {cap:.2f}'],
          dy=2.95, color=PLUCK)
    # scale bar
    ax.plot([L - 200, L], [1.2, 1.2], color=INK, lw=0.8)
    ax.text(L - 200 - 1.0 * upm, 1.2, '200 nt', ha='right', va='center', fontsize=FS)


# ------------------------------------------------------------------------------------------------ panel b
def lines(ax, x, y, rows, dy=3.4, ha='left', **kw):
    # one text per line on fixed baselines; y is the first baseline (mm)
    for i, r in enumerate(rows):
        ax.text(x, y - i * dy, r, ha=ha, va='baseline', fontsize=FS, **kw)


def box(ax, x, y, w, h, text, fc='white', ec=INK, lw=0.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=0.8', fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=FS, color=INK, linespacing=1.2)


def arrow(ax, x0, y0, x1, y1, color=INK):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='-|>', lw=0.6, color=color, mutation_scale=6,
                                                               shrinkA=0, shrinkB=0))


def panel_b(ax):
    # schematic, no data file; units are millimetres of the printed figure
    ax.set_xlim(0, W_MM); ax.set_ylim(0, ROW_B[1]); ax.axis('off')
    yt = ROW_B[1] - 1.0
    # (i) toy geometry: R = 12, k = 3, so each probe owns R - k + 1 = 10 read starts
    X0, R_, k_ = 8, 12, 3
    probes = [6, 16, 26, 46]; removed = 36
    ax.text(X0 - 2, yt, 'Each read contains at most one probe', ha='left', va='top', fontsize=7.5, fontweight='bold')
    yl = 43.0
    ax.plot([X0 - 3, X0 + 52], [yl, yl], color=INK, lw=0.8)
    for p in probes:
        ax.add_patch(Rectangle((X0 + p, yl - 0.9), k_, 1.8, fc=PLUCK, ec='white', lw=0.3, zorder=3))
    ax.add_patch(Rectangle((X0 + removed, yl - 0.9), k_, 1.8, fc='white', ec=REFU, lw=0.6, zorder=3))
    ax.text(X0 + removed + 1.5, yl + 1.5, 'removed (founder variant)', ha='center', va='bottom', fontsize=FS, color=REFU)
    for j, p in enumerate(probes):
        ax.add_patch(Rectangle((X0 + p - (R_ - k_), yl - 3.0), R_ - k_ + 1, 1.1, fc=PLUCK if j % 2 == 0 else GP, alpha=0.4, ec='none'))
    ax.annotate('', xy=(X0 + 7, yl - 4.6), xytext=(X0 + 16.95, yl - 4.6),
                arrowprops=dict(arrowstyle='<->', lw=0.5, color=INK, shrinkA=0, shrinkB=0, mutation_scale=4))
    ax.text(X0 + 17.8, yl - 4.6, '$R-k+1$ starts', ha='left', va='center', fontsize=FS)
    starts = [-2, 2, 5, 9, 13, 18, 21, 24, 28, 31, 34, 38, 41]
    for i, s0 in enumerate(starts):
        hit = [j for j, p in enumerate(probes) if p - (R_ - k_) <= s0 <= p]
        c = (PLUCK if hit[0] % 2 == 0 else GP) if hit else '#9E9E9E'
        y_ = yl - 7.4 - (i % 5) * 1.45
        ax.annotate('', xy=(X0 + s0 + R_, y_), xytext=(X0 + s0, y_),
                    arrowprops=dict(arrowstyle='-|>', lw=0.9, color=c, mutation_scale=4, shrinkA=0, shrinkB=0))
    lines(ax, X0 - 2, 25.0, ['Each read (arrow, length $R$) takes the colour', 'of the probe whose window holds its start;',
                             'grey reads start in a removed probe\'s window', 'and are not counted. $E[S_g]\\approx$ fragments of $g$',
                             'times $\\kappa_g = m_g(R-k+1)/\\tilde{l}_g$.'])
    # (ii) one streaming pass
    x2 = 66
    ax.text(x2, yt, 'One streaming pass per sample', ha='left', va='top', fontsize=7.5, fontweight='bold')
    box(ax, x2, 38.0, 17, 8.0, 'gzipped\nFASTQ')
    arrow(ax, x2 + 17, 42.0, x2 + 21, 42.0)
    box(ax, x2 + 21, 36.5, 33, 11.0, 'each 25-mer of each\nread looked up in\none probe table')
    arrow(ax, x2 + 37.5, 36.5, x2 + 37.5, 32.5)
    box(ax, x2 + 21, 26.0, 33, 6.5, 'probe counts $c_j$', ec=PLUCK)
    arrow(ax, x2 + 21, 30.6, x2 + 17.5, 30.6); arrow(ax, x2 + 21, 27.6, x2 + 17.5, 27.6)
    ax.text(x2 + 17, 30.6, '$S_g=\\Sigma_j c_j$', ha='right', va='center', fontsize=FS)
    ax.text(x2 + 17, 27.6, '$N_g^{\\mathrm{EM}}$ (EM k-mers)', ha='right', va='center', fontsize=FS)
    lines(ax, x2, 20.0, ['No alignment, no per-sample index: one', 'table per annotation, variant set and',
                         'read length serves every animal.'])
    # (iii) hybrid switch and design-only flag
    x3 = 124
    ax.text(x3, yt, 'Hybrid switch and reliability flag', ha='left', va='top', fontsize=7.5, fontweight='bold')
    cx, cy, hw, hh = x3 + 12, 41.5, 12.0, 6.5
    ax.add_patch(Polygon([[cx, cy + hh], [cx + hw, cy], [cx, cy - hh], [cx - hw, cy]], closed=True, fc='white', ec=INK, lw=0.6))
    ax.text(cx, cy, '$f_g \\geq \\tau$ ?', ha='center', va='center', fontsize=FS)
    box(ax, x3 + 26, 43.6, 28, 7.6, 'yes: gene probes\n$N_g = S_g/\\kappa_g$', ec=GP)
    box(ax, x3 + 26, 32.2, 28, 7.6, 'no: transcript EM over\ngene-unique k-mers', ec=EM)
    arrow(ax, cx + hw * 0.5, cy + hh * 0.5, x3 + 26, 47.4); arrow(ax, cx + hw * 0.5, cy - hh * 0.5, x3 + 26, 36.0)
    lines(ax, x3, 28.4, ['observed capture, per cohort:', '$f_g=\\Sigma_s S_{gs}/\\Sigma_s N^{\\mathrm{EM}}_{gs}$; $\\tau$ = 0.6'],
          dy=4.0)
    box(ax, x3, 5.0, 54, 15.0, '', fc=PALE, ec=MID)
    lines(ax, x3 + 27, 17.3, ['set from the design alone, before counting:', 'flag if $u_g$ < 0.6, where $u_g$ is the share',
                              'of EM-set windows that are gene-unique', 'and single-copy'], ha='center')


# ------------------------------------------------------------------------------------------------ panel c
CAPX = 'nominal capture $\\kappa_g = m_g(R-k+1)/\\tilde{l}_g$'


def panel_c(ax):
    # data: fig1c_design_capture.tsv
    d = read('fig1c_design_capture.tsv')
    styles = {'50 bp': dict(ls='-', lw=1.4), '101 bp': dict(ls=(0, (4, 1.5)), lw=1.4)}
    t, handles = {}, []
    for des, x in d.groupby('design', sort=False):
        q = np.sort(x.design_capture.to_numpy()); e = np.arange(1, len(q) + 1) / len(q)
        ax.step(q, e, where='post', color=PLUCK, **styles[des])
        t[des] = dict(n=int(x.gene_probes.sum()), low=(x.design_capture < 0.6).mean(), fl=int(x.flag.sum()), flp=x.flag.mean(),
                      zero=int((x.gene_probes == 0).sum()), sp=int(x.spacing_nt.iloc[0]))
        handles.append(Line2D([], [], color=PLUCK, **styles[des],
                              label=f'{des} design, {t[des]["sp"]} nt\n{t[des]["n"]:,} gene probes'))
    ax.axvline(0.6, color=MID, lw=0.6, ls=':')
    for des, dx in (('50 bp', -0.02), ('101 bp', 0.02)):
        ax.plot(0.6, t[des]['low'], 'o', ms=2.5, color=PLUCK, mec='white', mew=0.4, zorder=5)
        ax.text(0.6 + dx, t[des]['low'] + (0.01 if dx < 0 else -0.045), f'{t[des]["low"]:.0%}', ha='right' if dx < 0 else 'left',
                va='bottom' if dx < 0 else 'top', fontsize=FS, color=PLUCK)
    ax.text(0.615, 0.015, '0.6', ha='left', va='bottom', fontsize=FS, color=MID)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.0)
    ax.set_xlabel(CAPX); ax.set_ylabel('cumulative share of 16,885 core genes')
    a, b = t['50 bp'], t['101 bp']
    ax.text(0.03, 0.665, f'flagged ($u_g$ < 0.6):\n{a["flp"]:.1%} (50 bp)\n{b["flp"]:.1%} (101 bp)\nno gene probe: {a["zero"]}',
            transform=ax.transAxes, ha='left', va='top', fontsize=FS, linespacing=1.3)
    lg = ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.015, 1.02), frameon=True, handlelength=1.5,
                   handletextpad=0.5, labelspacing=0.5, borderpad=0.3)
    lg.get_frame().set(facecolor='white', edgecolor='none', alpha=1); lg._legend_box.align = 'left'


# ------------------------------------------------------------------------------------------------ panel d
YLIM = (5e-3, 20)


def panel_d(axs):
    # data: fig1d_model_check.tsv, fig1d_model_check_summary.tsv
    d = read('fig1d_model_check.tsv'); s = read('fig1d_model_check_summary.tsv').set_index('cohort')
    bins = np.logspace(np.log10(0.01), 0, 13)
    for ax, coh in zip(axs, ['Keele kidney, 50 bp', 'GSE199702 liver, 101 bp']):
        x = d[(d.cohort == coh) & d.used & (d.realized_capture > 0)]
        for br, c, mk, z in (('EM', EM, '^', 2), ('gene probes', GP, 'o', 3)):
            u = x[x.branch == br]
            # genes beyond the axis limits are drawn on the edge (kidney: 22 above, 36 below; liver: 45 below)
            ax.scatter(u.design_capture, np.clip(u.realized_capture, YLIM[0] * 1.03, YLIM[1] / 1.03), s=1.4, c=c, marker=mk, lw=0, alpha=0.4, zorder=z, rasterized=True)
        ax.plot([1e-3, 2], [1e-3, 2], color=INK, lw=0.7, zorder=4)
        # binned median and interquartile range over all used genes
        b = np.digitize(x.design_capture, bins)
        mids, med, lo, hi = [], [], [], []
        for i in range(1, len(bins)):
            v = x.realized_capture[b == i]
            if len(v) >= 20:
                mids.append(np.sqrt(bins[i - 1] * bins[i])); med.append(np.median(v)); lo.append(np.quantile(v, .25)); hi.append(np.quantile(v, .75))
        ax.errorbar(mids, med, yerr=[np.array(med) - lo, np.array(hi) - med], fmt='o', ms=2.6, color=INK, mfc='white', mew=0.7,
                    elinewidth=0.7, capsize=0, zorder=6)
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(3e-3, 1.4); ax.set_ylim(*YLIM)
        ax.set_title(coh, pad=3)
        ax.set_xlabel(CAPX)
        r = s.loc[coh]
        ax.text(0.03, 0.97, f'{int(r.genes):,} genes, {int(d[d.cohort == coh].n_samples.iloc[0])} samples\n'
                            f'median ratio {r.median_ratio:.2f}\nslope {r.theil_sen_slope:.2f} [{r.slope_lo:.2f}, {r.slope_hi:.2f}]',
                transform=ax.transAxes, ha='left', va='top', fontsize=FS, linespacing=1.3)
        ax.text(0.97, 0.03, f'gene probes {r["share_gene probes"]:.0%}\nEM {r.share_EM:.0%}',
                transform=ax.transAxes, ha='right', va='bottom', fontsize=FS, linespacing=1.3)
        ax.text(0.0115, 0.0068, '$y = x$', fontsize=FS, color=INK)
    axs[0].set_ylabel('realized capture $\\Sigma_s S_{gs}\\,/\\,\\Sigma_s\\,\\mathrm{truth}_{gs}$')
    axs[1].tick_params(labelleft=False)
    axs[1].legend(handles=[Line2D([], [], color=GP, marker='o', ls='none', ms=3, label='gene probes chosen'),
                           Line2D([], [], color=EM, marker='^', ls='none', ms=3, label='EM chosen'),
                           Line2D([], [], color=INK, marker='o', mfc='white', ls='none', ms=3, label='binned median, IQR')],
                  loc='upper left', bbox_to_anchor=(0.0, 0.72), frameon=False, handletextpad=0.2, labelspacing=0.3)


def main():
    fig = plt.figure(figsize=(W_MM * MM, H_MM * MM))
    f = lambda x, y, w, h: fig.add_axes([x / W_MM, y / H_MM, w / W_MM, h / H_MM])
    ax_a = f(0, ROW_A[0], W_MM, ROW_A[1])
    ax_b = f(0, ROW_B[0], W_MM, ROW_B[1])
    y, h = ROW_CD
    ax_c = f(13, y, 47, h)
    ax_d1 = f(80, y, 46, h); ax_d2 = fig.add_axes([131 / W_MM, y / H_MM, 46 / W_MM, h / H_MM], sharey=ax_d1)
    if os.path.exists(f'{D}/fig1a_gene_positions.tsv'):
        panel_a(ax_a, fig)
    panel_b(ax_b); panel_c(ax_c); panel_d([ax_d1, ax_d2])
    letter(fig, 'a', 1.0, H_MM - 1.0)
    letter(fig, 'b', 1.0, ROW_B[0] + ROW_B[1] + 1.5)
    letter(fig, 'c', 1.0, y + h + 6.5)
    letter(fig, 'd', 64.0, y + h + 6.5)
    for ext in ('pdf', 'svg'):
        fig.savefig(f'{OUT}/Fig1.{ext}', dpi=300)
    fig.savefig(PREVIEW, dpi=220)


if __name__ == '__main__':
    main()
