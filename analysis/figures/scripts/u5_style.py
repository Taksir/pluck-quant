"""Shared colours, markers and fonts for Fig. 5 and Figs S8, S9, S12."""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Folder constants below point to the original working folder; set them for your copy.
DATA = '/mnt/project-files/paper/figures/data'
FIGS = '/mnt/project-files/paper/figures'
MM = 1 / 25.4

PIPE = {  # label: (colour, marker, size factor, filled)
    'PLUCK': ('#D55E00', 'o', 1.4, True),
    'Salmon-pan': ('#0072B2', 's', 1.0, True),
    'Salmon-ref': ('#56B4E9', 's', 1.0, False),
    'STAR+fC': ('#009E73', '^', 1.0, True),
    # ablation arms: own filled shapes so no panel separates PLUCK and an arm by shade alone (as in fig2_style.py)
    'gene probes alone': ('#F08A5D', 'h', 1.0, True),
    'EM alone': ('#8C2D04', 'p', 1.0, True),
    'all gene-unique k-mers': ('#FDB863', 'P', 1.0, True),   # arm 2 (Fig. S12c)
    'flag-guided Salmon second pass': ('#000000', '*', 1.0, True),   # Fig. S12a, b
}
CODE2LABEL = {'hybrid': 'PLUCK', 'vf': 'PLUCK', 'salmon_pan': 'Salmon-pan', 'salmon_ref': 'Salmon-ref', 'star_fc': 'STAR+fC'}
COHORT_ORDER = ['keele_kidney', 'keele_liver', 'munger', 'gse199702_liver', 'gse199702_heart', 'gse199702_kidney']
COHORT_LABEL = {'keele_kidney': 'Keele kidney, 50 bp', 'keele_liver': 'Keele liver, 50 bp†', 'munger': 'Munger founders, 100 bp',
                'gse199702_liver': 'GSE199702 liver, 101 bp', 'gse199702_heart': 'GSE199702 heart, 101 bp', 'gse199702_kidney': 'GSE199702 kidney, 101 bp'}
COHORT_SHORT = {'keele_kidney': 'Keele\nkidney\n50 bp', 'keele_liver': 'Keele\nliver\n50 bp†', 'munger': 'Munger\nfounders\n100 bp',
                'gse199702_liver': 'GSE199702\nliver\n101 bp', 'gse199702_heart': 'GSE199702\nheart\n101 bp', 'gse199702_kidney': 'GSE199702\nkidney\n101 bp'}
GREY_B1 = '#9A9A9A'
GREY_B2 = '#4D4D4D'

def setup():
    # Elsevier artwork rule: lettering >= 7 pt and sub/superscripts >= 6 pt at final size. matplotlib shrinks mathtext scripts to 0.7 of the
    # base size (4.9 pt at 7 pt); 0.87 gives 6.1 pt. No rcParam exists for this, so the module constant is set.
    import matplotlib._mathtext as _mt
    _mt.SHRINK_FACTOR = 0.87
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
        'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'lines.linewidth': 0.9, 'axes.titlesize': 7, 'figure.dpi': 150})

def ax_mm(fig, x, y, w, h):
    """Axes placed in millimetres from the top-left corner of the figure."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])

def letter(fig, x, y, s):
    W, H = fig.get_size_inches() / MM
    fig.text(x / W, 1 - y / H, s, fontsize=10, fontweight='bold', va='top', ha='left')

def mtext(fig, x, y, s, **kw):
    W, H = fig.get_size_inches() / MM
    kw.setdefault('fontsize', 7)
    return fig.text(x / W, 1 - y / H, s, **kw)

EDGE = {'all gene-unique k-mers': '#B07020'}   # light fill needs a darker outline

def pmark(ax, x, y, label, size=4.0, z=3, face=None, scale=1.0, alpha=1.0, **kw):
    """One pipeline marker. A second, paired value of the same pipeline (other calibration, other gene set) keeps the
    pipeline's own marker and is drawn smaller and lighter (scale=0.75, alpha=0.55), as in Fig. 3c/d; fill never carries that meaning."""
    col, mk, f, filled = PIPE[label]
    fc = face if face is not None else (col if filled else 'white')
    ec = EDGE.get(label, 'white' if (label == 'PLUCK' and fc != 'white') else col)
    ax.plot(x, y, ls='none', marker=mk, ms=size * f * scale, mfc=fc, mec=ec, mew=0.6 if ec in ('white', '#B07020') else 0.9, alpha=alpha,
            zorder=z + (2 if label == 'PLUCK' else 0), **kw)

def save(fig, name):
    fig.savefig(f'{FIGS}/{name}.pdf')
    fig.savefig(f'{FIGS}/{name}.svg')
