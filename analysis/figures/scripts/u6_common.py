"""Shared colours, markers, line styles, labels and order for Fig. 6 and Figs S10-S11.

Salmon-ref 25% and 10% share the Salmon-ref colour and are told apart by line style plus a direct label.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
DATA = '/mnt/project-files/paper/figures/data'
OUT = '/mnt/project-files/paper/figures'

MM = 1 / 25.4
BASE_MS = 4.0  # marker size (points) for every pipeline except PLUCK

# label: (colour, marker, filled, size factor, linestyle, line width factor, order)
ENC = {
    'PLUCK':              ('#D55E00', 'o', True, 1.4, '-', 1.6, 1),
    'Salmon-pan':         ('#0072B2', 's', True, 1.0, '-', 1.0, 2),
    'Salmon-ref':         ('#56B4E9', 's', False, 1.0, '-', 1.0, 3),
    'STAR+fC':            ('#009E73', '^', True, 1.0, '-', 1.0, 4),
    'kallisto-pan':       ('#882255', 'D', True, 1.0, '-', 1.0, 5),
    'kallisto-ref':       ('#CC79A7', 'D', False, 1.0, '-', 1.0, 6),
    'Kmerator+countTags': ('#7F7F7F', 'v', True, 1.0, '-', 1.0, 7),
    'RNA-Skim':           ('#333333', 'x', True, 1.0, '-', 1.0, 8),
    'ref-unique probes':  ('#332288', 'o', False, 1.0, '--', 1.0, 9),
    'Salmon-pan 25%':     ('#4C94C8', 's', True, 0.8, '--', 1.0, 10),
    'Salmon-pan 10%':     ('#99C2DD', 's', True, 0.7, ':', 1.0, 11),
    'Salmon-ref 25%':     ('#56B4E9', 's', False, 0.8, '--', 1.0, 12),
    'Salmon-ref 10%':     ('#56B4E9', 's', False, 0.7, ':', 1.0, 13),
    'GBRS':               ('#B8B8B8', '^', False, 1.0, 'none', 1.0, 14),
}

# method codes in the source TSVs -> label
CODE = {
    'hybrid': 'PLUCK', 'salmon_pan': 'Salmon-pan', 'salmon_ref': 'Salmon-ref',
    'star_featurecounts': 'STAR+fC', 'star_fc': 'STAR+fC', 'kallisto_pan': 'kallisto-pan',
    'kallisto_ref': 'kallisto-ref', 'kmerator_counttags': 'Kmerator+countTags', 'kmerator': 'Kmerator+countTags',
    'rnaskim': 'RNA-Skim', 'ref_unique': 'ref-unique probes', 'ref': 'ref-unique probes',
    'salmon_pan_p25': 'Salmon-pan 25%', 'salmon_pan_p10': 'Salmon-pan 10%',
    'salmon_ref_p25': 'Salmon-ref 25%', 'salmon_ref_p10': 'Salmon-ref 10%',
    'pan_p25': 'Salmon-pan 25%', 'pan_p10': 'Salmon-pan 10%', 'ref_p25': 'Salmon-ref 25%', 'ref_p10': 'Salmon-ref 10%',
    'gbrs': 'GBRS', 'hybrid ': 'PLUCK',
}


def order(label):
    return ENC[label][6]


def mkw(label, ms_scale=1.0):
    """Marker keyword arguments for plot()/errorbar()."""
    c, m, filled, sf, ls, lwf, _ = ENC[label]
    kw = dict(marker=m, color=c, markersize=BASE_MS * sf * ms_scale, markeredgewidth=0.8)
    if label == 'PLUCK':
        kw.update(markerfacecolor=c, markeredgecolor='white', zorder=10)
    elif m == 'x':
        kw.update(markeredgecolor=c, markeredgewidth=1.1, zorder=4)
    elif filled:
        kw.update(markerfacecolor=c, markeredgecolor=c, zorder=4)
    else:
        kw.update(markerfacecolor='white', markeredgecolor=c, zorder=4)
    return kw


def lkw(label, base_lw=0.9):
    c, m, filled, sf, ls, lwf, _ = ENC[label]
    return dict(color=c, linestyle=('-' if ls == 'none' else ls), linewidth=base_lw * lwf,
                zorder=9 if label == 'PLUCK' else 3)


def read(name):
    return pd.read_csv(f'{DATA}/{name}', sep='\t', comment='#')


def setup():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'xtick.labelsize': 7,
        'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.linewidth': 0.6,
        'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4,
        'ytick.minor.width': 0.4, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
        'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'axes.titlepad': 3, 'axes.labelpad': 2,
    })


def panel_letter(fig, ax, letter, dx=-0.0, dy=0.0):
    """Bold lowercase letter at the top-left corner of an axes (figure coordinates)."""
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=10, fontweight='bold', va='bottom', ha='left')


def save(fig, stem):
    fig.savefig(f'{OUT}/{stem}.pdf')
    fig.savefig(f'{OUT}/{stem}.svg')
