"""Shared colours, markers and fonts for Fig. 2 and Figs S1-S3."""
import matplotlib as mpl
import matplotlib._mathtext as _mathtext

# Folder constants below point to the original working folder; set them for your copy.
DATA = '/mnt/project-files/paper/figures/data'
FIGS = '/mnt/project-files/paper/figures'

# label, colour, marker, filled, size factor, line style
P = {
    'hybrid':      ('PLUCK', '#D55E00', 'o', True, 1.4, '-'),
    'salmon_pan':  ('Salmon-pan', '#0072B2', 's', True, 1.0, '-'),
    'salmon_ref':  ('Salmon-ref', '#56B4E9', 's', False, 1.0, '-'),
    'star_fc':     ('STAR+fC', '#009E73', '^', True, 1.0, '-'),
    'kallisto_pan': ('kallisto-pan', '#882255', 'D', True, 1.0, '-'),
    'kallisto_ref': ('kallisto-ref', '#CC79A7', 'D', False, 1.0, '-'),
    'kmerator':    ('Kmerator+countTags', '#7F7F7F', 'v', True, 1.0, '-'),
    'rnaskim':     ('RNA-Skim', '#333333', 'x', True, 1.0, '-'),
    'ref':         ('ref-unique probes', '#332288', 'o', False, 1.0, '--'),
    # ablation arms
    # distinct shapes so no panel separates PLUCK and its branches by colour alone
    'vf':          ('gene probes alone', '#F08A5D', 'h', True, 1.0, '-'),
    'emg':         ('EM alone', '#8C2D04', 'p', True, 1.0, '-'),
}
ORDER = ['hybrid', 'salmon_pan', 'salmon_ref', 'star_fc', 'kallisto_pan', 'kallisto_ref', 'kmerator', 'rnaskim', 'ref']
LAYOUT = {'se50': 'Simulation, 50 bp SE', 'se101': 'Simulation, 101 bp SE', 'pe': 'Simulation, 101 bp PE'}
FOUNDER = {'A': 'A/J', 'C': '129S1', 'D': 'NOD', 'E': 'NZO', 'F': 'CAST', 'G': 'PWK', 'H': 'WSB'}   # standard CC order, B6 omitted
FOUNDERS = list('ACDEFGH')
GREY = '#666666'
FS_ANN = 7     # smallest text in these figures; Elsevier minimum 7 pt at final size
MM = 1 / 25.4


def setup():
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'legend.fontsize': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'axes.titlepad': 3,
        'savefig.dpi': 300, 'lines.linewidth': 0.9,
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic',
        'mathtext.bf': 'Liberation Sans:bold'})
    # sub/superscripts at 0.86 of the base size, so 7 pt text carries 6.0 pt scripts (Elsevier minimum 6 pt)
    _mathtext.SHRINK_FACTOR = 0.86


def mstyle(code, base=14):
    """scatter kwargs for a pipeline code; PLUCK gets the white edge and the largest marker."""
    lab, col, mk, filled, fac, _ = P[code]
    kw = dict(marker=mk, s=base * fac ** 2, linewidths=0.8, zorder=5 if code == 'hybrid' else 3)   # size factor applies to the diameter
    if mk == 'x':
        kw.update(c=col)
    elif filled:
        kw.update(facecolors=col, edgecolors='white' if code == 'hybrid' else col)
        if code == 'hybrid':
            kw.update(linewidths=0.35)
    else:
        kw.update(facecolors='white', edgecolors=col)
    return kw


def panel_letter(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')


def save(fig, name):
    for ext in ('pdf', 'svg'):
        fig.savefig(f'{FIGS}/{name}.{ext}', bbox_inches='tight', pad_inches=0.02)
