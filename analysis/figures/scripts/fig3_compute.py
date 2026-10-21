#!/usr/bin/env python3
"""Per-gene and per-sample values behind Fig. 3a, 3b, 3e and Fig. S4 (scoring set with all methods, real cohorts).

Reads (read-only) the scoring inputs and tool outputs of the comparison, through the gather script's own
loading functions, and replays the scorer's sample-gene filter and calibration:
  revision/final_comparison/frames/<cohort>/            (records, universe, samples; frames.sha256)
  revision/final_comparison/results_ll/acc/<cohort>/<srr>/*.genes.tsv.gz   (kallisto, Kmerator+countTags, STAR+fC)
  scoring/fc_gather_score.py and scoring/score_e9.py of this repository (imported; the scorer hash is checked on import)
Writes to figures/data/:
  fig3ab_pergene_u3.tsv        per gene: within-gene Spearman (rho_<tool>) and within-gene noise SD (sd_<tool>), flag
  fig3e_persample_u3.tsv       per sample: Spearman across genes against the cohort truth, and its rank among tools
  figS4_pairwise_u3.tsv        per sample and tool pair: Spearman across genes between the two estimates
  fig3_reproduction_u3.tsv     medians of the per-gene values against the scorer's numbers for the set with all methods
The per-gene definitions are the scorer's: rho = Spearman over the gene's sample-genes (>= 10) of log2(C*est+1) against
log2(truth+1); noise SD = 1.4826 * MAD over samples of log2((C*est+1)/(truth+1)) (>= 10 sample-genes). Sample-genes:
truth >= 50, stratum B6 or nonB6, an estimate from every tool of the set, calibration constant present (B6 rule in Keele,
non-hit rule in Munger).
Run: python3 -I figures/scripts/fig3_compute.py
"""
import collections, importlib.util, os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision'
FC = f'{REV}/final_comparison'
OUT = '/mnt/project-files/paper/figures/data'
FRAMES = f'{FC}/frames'
TOOLS_ROOT = f'{FC}/results_ll/acc'
LABEL = {'vf': 'hybrid', 'salmon_pan': 'Salmon-pan', 'salmon_ref': 'Salmon-ref', 'star_fc': 'STAR+fC',
         'kallisto_pan': 'kallisto-pan', 'kallisto_ref': 'kallisto-ref', 'kmerator': 'Kmerator+countTags',
         'ref': 'ref-unique probes'}
ORDER = ['vf', 'salmon_pan', 'salmon_ref', 'star_fc', 'kallisto_pan', 'kallisto_ref', 'kmerator', 'ref']

spec = importlib.util.spec_from_file_location('fcg', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scoring', 'fc_gather_score.py'))
fcg = importlib.util.module_from_spec(spec); spec.loader.exec_module(fcg)
fcg.load_scorer()                       # exits if the scorer's hash differs
S9 = fcg.S9


def run(name):
    c = fcg.load_cohort(name, FRAMES)
    acc = fcg.acc_for(name)
    fcg.attach_tools(c, TOOLS_ROOT, acc)
    fcg.set_calibration(name)
    methods = [m for m in ORDER if m in fcg.EXISTING + acc]
    cov = fcg.covered_samples(c, methods)
    subsets, _ = fcg.flag_sets(c)
    C = S9.calibrate(c.recs, methods)
    per_gene = collections.defaultdict(list)       # gene -> [(srr, lt, le[], lr[])]
    per_samp = collections.defaultdict(list)       # srr -> [(lt, le[])]
    for (srr, g), r in c.recs.items():
        if srr not in cov or r['truth'] < S9.MIN_TRUTH or r['stratum'] not in ('B6', 'nonB6'):
            continue
        if any(r['est'].get(m) is None or (srr, m) not in C for m in methods):
            continue
        le = [float(np.log2(C[(srr, m)] * r['est'][m] + 1)) for m in methods]
        lr = [float(np.log2((C[(srr, m)] * r['est'][m] + 1) / (r['truth'] + 1))) for m in methods]
        lt = float(np.log2(r['truth'] + 1))
        per_gene[g].append((lt, le, lr)); per_samp[srr].append((lt, le))
    # per gene
    rows = []
    for g, ps in per_gene.items():
        if len(ps) < S9.MIN_SPEARMAN_N:
            continue
        lt = np.array([p[0] for p in ps]); LE = np.array([p[1] for p in ps]); LR = np.array([p[2] for p in ps])
        row = dict(cohort=name, gene=g, n_samples=len(ps),
                   flag='flagged' if g in subsets['flagged'] else ('unflagged' if g in subsets['unflagged'] else 'not_in_flag_table'))
        for j, m in enumerate(methods):
            le = LE[:, j]
            row['rho_' + LABEL[m]] = float(spearmanr(le, lt).statistic) if (lt.std() > 0 and le.std() > 0) else np.nan
            x = LR[:, j]
            row['sd_' + LABEL[m]] = 1.4826 * float(np.median(np.abs(x - np.median(x))))
        rows.append(row)
    pg = pd.DataFrame(rows)
    # per sample: Spearman across genes against truth; rank among tools (1 = highest); pairwise between tools
    srows, prow = [], []
    for s, ps in per_samp.items():
        lt = np.array([p[0] for p in ps]); LE = np.array([p[1] for p in ps])
        M = np.column_stack([lt, LE])
        R = spearmanr(M).statistic
        rho = R[0, 1:]
        rk = rankdata(-rho, method='min')
        for j, m in enumerate(methods):
            srows.append(dict(cohort=name, srr=s, n_genes=len(ps), tool=LABEL[m], spearman_vs_truth=float(rho[j]),
                              rank=int(rk[j]), n_tools=len(methods)))
        names = ['cohort truth'] + [LABEL[m] for m in methods]
        for i in range(len(names)):
            for k in range(i + 1, len(names)):
                prow.append(dict(cohort=name, srr=s, tool_a=names[i], tool_b=names[k], spearman=float(R[i, k])))
    return pg, pd.DataFrame(srows), pd.DataFrame(prow), methods, len(cov)


def main():
    head = ('# source: revision/final_comparison/frames/<cohort>/ and results_ll/acc/ via fc_gather_score.py (70bdc836...a9f975) '
            'and score_e9.py (7710d525...467d); tools of the scoring set with all methods; scorer sample-gene filter and calibration replayed\n'
            '# written by: fig3_compute.py; values computed, not copied\n'
            '# note: truths are the cohort truths (RSEM on strain genomes), which share reads with full-read methods. '
            '"hybrid" is plotted as PLUCK.\n')
    pgs, sps, pps, rep = [], [], [], []
    ref = pd.read_csv(f'{OUT}/fig3_metrics_all_frames.tsv', sep='\t', comment='#')
    ref = ref[(ref.frame == 'A') & (ref.subset == 'all')]
    for name in fcg.COHORTS:
        pg, sp, pp, methods, ns = run(name)
        pgs.append(pg); sps.append(sp); pps.append(pp)
        for m in methods:
            t = LABEL[m]; r = ref[(ref.cohort == name) & (ref.tool == t)].iloc[0]
            rho = pg['rho_' + t].dropna(); sd = pg['sd_' + t]
            rep.append(dict(cohort=name, tool=t, n_samples=ns, median_rho=round(float(rho.median()), 4), n_rho=len(rho),
                            scorer_rho=r.within_gene_spearman, scorer_n_rho=r.n_genes_spearman,
                            median_sd=round(float(sd.median()), 4), n_sd=len(sd), scorer_sd=r.noise_sd, scorer_n_sd=r.n_genes_noise,
                            match=bool(round(float(rho.median()), 4) == r.within_gene_spearman and round(float(sd.median()), 4) == r.noise_sd)))
        print(name, 'done', flush=True)
    for fn, df, note in (('fig3ab_pergene_u3.tsv', pd.concat(pgs), '# file: fig3ab_pergene_u3.tsv (one row per gene with >= 10 sample-genes in the set with all methods)\n'),
                         ('fig3e_persample_u3.tsv', pd.concat(sps), '# file: fig3e_persample_u3.tsv (per-sample Spearman across genes, a different metric from within-gene Spearman)\n'),
                         ('figS4_pairwise_u3.tsv', pd.concat(pps), '# file: figS4_pairwise_u3.tsv (per-sample Spearman across genes between two estimates)\n'),
                         ('fig3_reproduction_u3.tsv', pd.DataFrame(rep), '# file: fig3_reproduction_u3.tsv (median of the per-gene values against the scorer, set with all methods, all genes)\n')):
        with open(f'{OUT}/{fn}', 'w') as f:
            f.write(note + head)
            df.to_csv(f, sep='\t', index=False, float_format='%.6g')
    print(pd.DataFrame(rep).to_string())


if __name__ == '__main__':
    main()
