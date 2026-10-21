"""Shared colours, markers and fonts for Fig. 3 and Figs S4-S7.
Data files are read only from figures/data/. Every panel function names its data file in a comment."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
DATA = '/mnt/project-files/paper/figures/data'
FIGS = '/mnt/project-files/paper/figures'
MM = 1 / 25.4

# label, colour, marker, filled, size factor (order = list order)
PIPE = [
    ('PLUCK', '#D55E00', 'o', True, 1.4),
    ('Salmon-pan', '#0072B2', 's', True, 1.0),
    ('Salmon-ref', '#56B4E9', 's', False, 1.0),
    ('STAR+fC', '#009E73', '^', True, 1.0),
    ('kallisto-pan', '#882255', 'D', True, 0.9),
    ('kallisto-ref', '#CC79A7', 'D', False, 0.9),
    ('Kmerator+countTags', '#7F7F7F', 'v', True, 1.0),
    ('RNA-Skim', '#333333', 'x', True, 1.0),
    ('ref-unique probes', '#332288', 'o', False, 1.0),
]
ORDER = [p[0] for p in PIPE]
STYLE = {p[0]: dict(color=p[1], marker=p[2], filled=p[3], size=p[4]) for p in PIPE}
TSV2LABEL = {'hybrid': 'PLUCK'}            # TSV tool value -> figure label; others unchanged

# cohort labels
COHORT_LABEL = {'keele_kidney': 'Keele kidney, 50 bp', 'keele_liver': 'Keele liver, 50 bp†',
                'munger': 'Munger founders, 100 bp'}
COHORTS = ['keele_kidney', 'keele_liver', 'munger']
NAME2KEY = {'Keele kidney': 'keele_kidney', 'Keele liver': 'keele_liver', 'Munger founders': 'munger'}


def setup():
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'legend.fontsize': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
        'xtick.major.size': 2, 'ytick.major.size': 0, 'xtick.minor.size': 1.2, 'xtick.minor.width': 0.4,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'axes.spines.top': False,
        'axes.spines.right': False, 'lines.linewidth': 0.8, 'axes.titlepad': 2, 'xtick.major.pad': 1.5,
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic',
        'mathtext.bf': 'Liberation Sans:bold'})


def read(name):
    return pd.read_csv(f'{DATA}/{name}', sep='\t', comment='#')


def label(tool):
    return TSV2LABEL.get(tool, tool)


def mk(ax, x, y, tool, scale=1.0, alpha=1.0, z=None, ms=3.2):
    """Plot one point in the pipeline's encoding. PLUCK is drawn on top with a white edge."""
    s = STYLE[tool]
    size = ms * s['size'] * scale
    kw = dict(marker=s['marker'], ms=size, alpha=alpha, ls='none', zorder=z or (5 if tool == 'PLUCK' else 3))
    if s['marker'] == 'x':
        kw.update(color=s['color'], mew=0.9)
    elif s['filled']:
        kw.update(mfc=s['color'], mec='white' if tool == 'PLUCK' else s['color'], mew=0.5 if tool == 'PLUCK' else 0.4)
    else:
        kw.update(mfc='white', mec=s['color'], mew=0.9)
    return ax.plot([x], [y], **kw)


def panel_letter(fig, x_mm, y_mm, letter, H_mm):
    fig.text(x_mm * MM / fig.get_figwidth(), 1 - y_mm * MM / fig.get_figheight(), letter, fontsize=10,
             fontweight='bold', va='top', ha='left')


def ax_mm(fig, x, y, w, h, W, H):
    """Axes from millimetres measured from the top-left corner."""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def save(fig, stem):
    fig.savefig(f'{FIGS}/{stem}.pdf')
    fig.savefig(f'{FIGS}/{stem}.svg')
    print('wrote', f'{FIGS}/{stem}.pdf', f'{FIGS}/{stem}.svg')
