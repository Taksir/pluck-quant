"""Shared helpers: harmonised design features per cohort layout, metrics."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
LAYOUT_FEAT = {'50SE': ('feat26.tsv.gz', '26'), '100SE': ('feat76.tsv.gz', '76'), '101SE': ('feat76.tsv.gz', '77')}
BASE = ['ntx', 'ntx_pc', 'ntx_ri_nmd', 'ntsl1', 'tsl1_len_logratio', 'ref_len', 'nB', 'frac_shared', 'frac_multicopy', 'frac_emg',
        'n_emg', 'frac_tx_invisible', 'frac_pc_invisible', 'em_capture', 'em_nvis', 'em_rank_def', 'em_var_factor', 'em_ambiguity',
        'em_npat', 'em_null_frac', 'em_log10_cond']
SPACED = ['m_vf', 'm_ref', 'vf_loss', 'gp_capture', 'f_pred_all', 'f_pred_pc']
LOG = {'ntx', 'ntx_pc', 'ntx_ri_nmd', 'ntsl1', 'ref_len', 'nB', 'n_emg', 'em_nvis', 'em_npat', 'm_vf', 'm_ref', 'em_rank_def'}
def design_features(featdir, layout):
    fn, sp = LAYOUT_FEAT[layout]
    F = pd.read_csv(f'{featdir}/{fn}', sep='\t')
    X = F[['gene_id'] + BASE].copy()
    for c in SPACED: X[c] = F[f'{c}_{sp}']
    for c in LOG: X[c] = np.log1p(X[c])
    X['em_var_factor'] = np.log2(X.em_var_factor.clip(lower=0.1, upper=1000))
    X['em_ambiguity'] = np.log1p(X.em_ambiguity.clip(lower=0, upper=100))
    for c in ('f_pred_all', 'f_pred_pc'): X[c] = X[c].clip(0, 3)
    X['pred_em_all'] = (F[f'f_pred_all_{sp}'] < 0.6) | F[f'f_pred_all_{sp}'].isna() | (F[f'm_vf_{sp}'] == 0)
    X['pred_em_pc'] = (F[f'f_pred_pc_{sp}'] < 0.6) | F[f'f_pred_pc_{sp}'].isna() | (F[f'm_vf_{sp}'] == 0)
    X['no_gene_probe'] = F[f'm_vf_{sp}'] == 0
    for c in ('pred_em_all', 'pred_em_pc', 'no_gene_probe'): X[c] = X[c].astype(float)
    return X
FEATS = BASE + SPACED + ['pred_em_all', 'pred_em_pc', 'no_gene_probe']
def auc(y, s):
    m = np.isfinite(s); y = np.asarray(y)[m]; s = np.asarray(s)[m]
    return (roc_auc_score(y, s) if 0 < y.sum() < len(y) else np.nan), (average_precision_score(y, s) if y.sum() else np.nan)
def at_flag(y, flag):
    y = np.asarray(y, bool); f = np.asarray(flag, bool)
    tp = (y & f).sum()
    return dict(flag_rate=f.mean(), precision=tp / max(f.sum(), 1), recall=tp / max(y.sum(), 1),
                tail_flagged=y[f].mean() if f.any() else np.nan, tail_unflagged=y[~f].mean() if (~f).any() else np.nan,
                enrichment=(y[f].mean() / y[~f].mean()) if f.any() and y[~f].mean() > 0 else np.nan)
