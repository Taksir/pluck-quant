#!/usr/bin/env python3
"""Source data for Fig. 1. Writes TSVs into figures/data/; fig1.py draws from those files only.

Panel a  fig1a_gene_positions.tsv, fig1a_gene_variants.tsv, fig1a_gene_isoforms.tsv
         One real core gene, Slc25a47 (ENSMUSG00000048856, chr12, + strand, 2 TSL1 protein-coding isoforms), with the
         101 bp design (R - k + 1 = 77). The candidate pool is not stored, so it is rebuilt here
         with the rules of e9rs.py prep() (25-mers of the longest TSL1 isoform; present in every
         TSL1 isoform; >= 100 nt from both ends; canonical k-mer in no other gene of the Ensembl 102 cDNA + ncRNA set;
         <= 1 copy in GRCm38; founder mask = overlap with any MGP v5 SNP or indel call that is not 0/0 in that founder;
         greedy spacing >= 77 nt). The rebuilt ref-unique and PLUCK probe positions and masks must equal the stored
         file revision/respacing/probes/probes_s77.npz for this gene, or the script stops.
         Inputs: Ensembl release 102 GTF, cDNA, ncRNA and GRCm38 primary assembly (ftp.ensembl.org), MGP REL-1505
         strain VCFs read by tabix region queries (ftp.ebi.ac.uk). Pass --ens DIR with the four downloaded files.
Panel c  fig1c_design_capture.tsv: design capture q_g = m_g (R - k + 1) / l~_g for every core-tier gene at the 50 bp
         (26 nt spacing) and 101 bp (77 nt spacing) designs; m_g = PLUCK gene probes (in_vf), l~_g = mean over TSL1
         isoforms of max(L - R + 1, 25) (as in the PLUCK estimator); flag u_g < 0.6 from
         reliability_flag/out/flag_table.tsv.
Panel d  fig1d_model_check.tsv: per gene, pooled observed capture sum_s S_gs / sum_s truth_gs against q_g.
         Keele kidney 50 bp (held out): S from reliability_flag/out/core_keele_kidney_rebuilt.npz (stored gene-probe
         sums from e4_kidney/counts), truth from final_comparison/frames/keele_kidney/records.tsv.gz (GEO RSEM).
         GSE199702 liver 101 bp (held out, 53 scored samples): S = vf26 * m_g * 77 / l~_g, the PLUCK gene-probe
         estimator inverted exactly (e4b_full/tab/liver/core_all.tsv.gz), truth = GBRS/EMASE column of that table.
         Genes: >= 1 gene probe and median truth >= 50 (the scorer's expression floor). Branch: gene-probe branch if
         the hybrid used the gene-probe value (fo >= 0.6 in kidney; hyb == vf26 in GSE199702 liver, read on the samples
         where vf26 and emg differ).
usage: fig1_data.py [--ens DIR] [--skip-a]
"""
import argparse, gzip, json, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision'
OUT = '/mnt/project-files/paper/figures/data'
K = 25
GENES = json.load(open(f'{REV}/fast_competitors/e9_prep_cache/ref/genes.json'))
GID = np.array([g['gene_id'] for g in GENES]); GNAME = np.array([g['gene_name'] for g in GENES]); NG = len(GENES)
FOUNDERS = [('A', 'A_J'), ('C', '129S1_SvImJ'), ('D', 'NOD_ShiLtJ'), ('E', 'NZO_HlLtJ'),
            ('F', 'CAST_EiJ'), ('G', 'PWK_PhJ'), ('H', 'WSB_EiJ')]  # bit order of the stored founder mask
MGP = 'https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs'


def header(path, source, notes):
    with open(path, 'w') as f:
        f.write(f'# file: {os.path.basename(path)}\n# source: {source}\n'
                f'# written by: fig1_data.py\n')
        for n in notes:
            f.write(f'# note: {n}\n')


def ltilde(R):
    return np.array([np.mean([max(L - R + 1, K) for L in g['tx_lens']]) for g in GENES])


def gene_probe_counts(fname):
    P = dict(np.load(f'{REV}/respacing/probes/{fname}'))
    return np.bincount(P['gene'][P['in_vf']], minlength=NG), P


# ---------------------------------------------------------------- panel c
def panel_c():
    fl = pd.read_csv(f'{REV}/reliability_flag/out/flag_table.tsv', sep='\t')
    tiers = pd.read_csv(f'{REV}/frozen_v2/gene_tiers.tsv', sep='\t')
    core = set(tiers.loc[tiers.tier == 'core', 'gene_id'])
    assert fl.gene_id.isin(core).all() and len(fl) == len(core) == 16885
    gpos = {g: i for i, g in enumerate(GID)}
    gi = fl.gene_id.map(gpos).astype(int).to_numpy()
    rows = []
    for R, fn, flagcol in ((50, 'probes_s26.npz', 'flag_50bp'), (101, 'probes_s77.npz', 'flag_100bp')):
        m, _ = gene_probe_counts(fn)
        q = m[gi] * (R - K + 1) / ltilde(R)[gi]
        rows.append(pd.DataFrame({'design': f'{R} bp', 'spacing_nt': R - K + 1, 'gene_id': fl.gene_id,
                                  'gene_name': fl.gene_name, 'gene_probes': m[gi], 'ltilde': ltilde(R)[gi].round(2),
                                  'design_capture': q.round(5), 'flag': fl[flagcol].astype(int)}))
    d = pd.concat(rows)
    path = f'{OUT}/fig1c_design_capture.tsv'
    header(path, 'respacing/probes/probes_s26.npz, probes_s77.npz (in_vf); fast_competitors/e9_prep_cache/ref/genes.json '
                 '(TSL1 isoform lengths); reliability_flag/out/flag_table.tsv (flag_50bp, flag_100bp)',
           ['derived: design_capture = gene_probes * spacing / ltilde; ltilde = mean over TSL1 isoforms of max(L - R + 1, 25)',
            'the 100 bp flag (76 nt design) is used for the 101 bp row, as in flag_table.tsv'])
    d.to_csv(path, sep='\t', index=False, mode='a')
    for des, x in d.groupby('design', sort=False):
        print('c', des, 'genes', len(x), 'probes', int(x.gene_probes.sum()), 'median m', np.median(x.gene_probes),
              'zero', int((x.gene_probes == 0).sum()), 'q<0.6', round(float((x.design_capture < 0.6).mean()), 4),
              'median q', round(float(np.median(x.design_capture)), 3), 'flagged', int(x.flag.sum()))


# ---------------------------------------------------------------- panel d
def summarize(tag, S, T, q, branch):
    ok = (q > 0) & (np.nanmedian(T, 1) >= 50) & np.isfinite(T).all(1) & np.isfinite(S).all(1)
    obs = S.sum(1) / T.sum(1)
    return ok, obs


def panel_d():
    rows, stats_rows = [], []
    # Keele kidney, 50 bp
    R = 50; W = R - K + 1; lt = ltilde(R); m, _ = gene_probe_counts('probes_s26.npz')
    E = np.load(f'{REV}/reliability_flag/out/core_keele_kidney_rebuilt.npz')
    srrs = list(E['srrs']); S_all = E['S']; fo = E['fo']
    # the stored vf26 estimate must equal the PLUCK gene-probe estimator applied to the stored sums
    with np.errstate(invalid='ignore', divide='ignore'):
        chk = np.where(m > 0, lt / (np.maximum(m, 1) * W) * S_all, np.nan)
    assert np.allclose(chk, E['vf26'], equal_nan=True, rtol=1e-9), 'vf26 != S * l~ / (m W)'
    rec = pd.read_csv(f'{REV}/final_comparison/frames/keele_kidney/records.tsv.gz', sep='\t', usecols=['srr', 'gene', 'truth'])
    tr = rec.pivot_table(index='gene', columns='srr', values='truth')
    cnt = pd.Series(GNAME).value_counts(); uniq = set(cnt[cnt == 1].index)
    gi_of = {n: i for i, n in enumerate(GNAME) if n in uniq}
    tr = tr[[g in gi_of for g in tr.index]]
    gi = np.array([gi_of[g] for g in tr.index]); si = [srrs.index(s) for s in tr.columns]
    S = S_all[np.ix_(si, gi)].T; T = tr.to_numpy(); q = m[gi] * W / lt[gi]
    ok, obs = summarize('kidney', S, T, q, None)
    rows.append(pd.DataFrame({'cohort': 'Keele kidney, 50 bp', 'gene_id': GID[gi], 'gene_name': GNAME[gi], 'n_samples': T.shape[1],
                              'gene_probes': m[gi], 'design_capture': q, 'probe_sum': S.sum(1), 'truth_sum': T.sum(1),
                              'realized_capture': obs, 'branch': np.where(fo[gi] >= 0.6, 'gene probes', 'EM'),
                              'used': ok}))
    # GSE199702 liver, 101 bp
    R = 101; W = R - K + 1; lt = ltilde(R); m, _ = gene_probe_counts('probes_s77.npz')
    t = pd.read_csv(f'{REV}/e4b_full/tab/liver/core_all.tsv.gz', sep='\t', usecols=['srr', 'gene', 'truth', 'vf26', 'hyb', 'emg'])
    gpos = {g: i for i, g in enumerate(GID)}
    t['gi'] = t.gene.map(gpos); t = t.dropna(subset=['gi']); t['gi'] = t.gi.astype(int)
    t['S'] = t.vf26 * m[t.gi] * W / lt[t.gi]
    # branch read off the stored hybrid on the samples where the two branch values differ (both are 0 in some samples)
    t['gp'] = np.where(np.isclose(t.vf26, t.emg, rtol=1e-12, atol=0), np.nan,
                       (np.isclose(t.hyb, t.vf26, rtol=1e-12, atol=0) & t.vf26.notna()).astype(float))
    Sw = t.pivot_table(index='gi', columns='srr', values='S'); Tw = t.pivot_table(index='gi', columns='srr', values='truth')
    Tw = Tw.reindex(index=Sw.index, columns=Sw.columns)
    gi = Sw.index.to_numpy(); S = Sw.to_numpy(); T = Tw.to_numpy(); q = m[gi] * W / lt[gi]
    gp = t.groupby('gi').gp.mean().reindex(gi).to_numpy()
    ok, obs = summarize('liver101', S, T, q, None)
    rows.append(pd.DataFrame({'cohort': 'GSE199702 liver, 101 bp', 'gene_id': GID[gi], 'gene_name': GNAME[gi], 'n_samples': T.shape[1],
                              'gene_probes': m[gi], 'design_capture': q, 'probe_sum': S.sum(1), 'truth_sum': T.sum(1),
                              'realized_capture': obs, 'branch': np.where(gp == 1, 'gene probes', np.where(gp == 0, 'EM', 'mixed')),
                              'used': ok}))
    d = pd.concat(rows)
    for c in ('design_capture', 'realized_capture'):
        d[c] = d[c].round(6)
    d['probe_sum'] = d.probe_sum.round(3)
    for coh, x in d.groupby('cohort', sort=False):
        u = x[x.used & (x.realized_capture > 0)]
        ratio = u.realized_capture / u.design_capture
        med = float(np.median(ratio))
        ts = stats.theilslopes(np.log10(u.realized_capture), np.log10(u.design_capture))
        rho = stats.spearmanr(u.design_capture, u.realized_capture)[0]
        row = dict(cohort=coh, genes=int(x.used.sum()), genes_zero_obs=int((x.used & (x.realized_capture <= 0)).sum()),
                   median_ratio=round(med, 4), theil_sen_slope=round(ts.slope, 4), slope_lo=round(ts.low_slope, 4),
                   slope_hi=round(ts.high_slope, 4), spearman=round(rho, 4))
        for b in ('gene probes', 'EM'):
            ub = u[u.branch == b]; r0 = np.log2(ub.realized_capture / ub.design_capture)
            row[f'n_{b}'] = len(ub)
            row[f'share_{b}'] = round(len(ub) / len(u), 4)
            row[f'median_q_{b}'] = round(float(np.median(ub.design_capture)), 3)
            row[f'iqr_log2_ratio_{b}'] = f'{np.quantile(r0, .25):.3f} to {np.quantile(r0, .75):.3f}'
            row[f'iqr_log2_ratio_over_median_{b}'] = f'{np.quantile(r0 - np.log2(med), .25):.3f} to {np.quantile(r0 - np.log2(med), .75):.3f}'
        row['mixed'] = int((u.branch == 'mixed').sum())
        stats_rows.append(row); print('d', row)
    path = f'{OUT}/fig1d_model_check.tsv'
    header(path, 'reliability_flag/out/core_keele_kidney_rebuilt.npz (S, fo); final_comparison/frames/keele_kidney/records.tsv.gz '
                 '(truth); e4b_full/tab/liver/core_all.tsv.gz (vf26, hyb, truth); respacing/probes/probes_s26.npz, probes_s77.npz',
           ['derived: realized_capture = probe_sum / truth_sum, summed over the cohort samples (not the switch statistic f_g, which uses N^EM); '
             'design_capture = nominal capture kappa_g = m (R - k + 1) / l~ (Methods N2.1)',
            'GSE199702 liver probe_sum is the PLUCK gene-probe estimator inverted: vf26 * m * 77 / l~ (exact)',
            'used = at least one gene probe, all samples with truth, median truth >= 50'])
    d.to_csv(path, sep='\t', index=False, mode='a')
    path = f'{OUT}/fig1d_model_check_summary.tsv'
    header(path, 'fig1d_model_check.tsv', ['derived: median_ratio = median over used genes of realized/nominal capture; '
           'theil_sen_slope of log10 realized on log10 nominal capture with its 95% interval; '
           'iqr_log2_ratio = IQR of log2(realized/nominal); iqr_log2_ratio_over_median = the same divided by median_ratio',
           'genes with zero probe sum are excluded from the slope and the ratio (count in genes_zero_obs)'])
    pd.DataFrame(stats_rows).to_csv(path, sep='\t', index=False, mode='a')


# ---------------------------------------------------------------- panel a
LUT = np.full(256, 4, np.uint8)
for i, ch in enumerate(b'ACGT'):
    LUT[ch] = i; LUT[ch + 32] = i


def codes(b):
    """2-bit forward and reverse-complement codes of every k-mer; valid = no N (as e9rs.seqcodes)."""
    n = len(b) - K + 1
    if n <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, np.uint64), np.zeros(0, bool)
    bb = b.astype(np.uint64); bad = (b == 4)
    f = np.zeros(n, np.uint64); r = np.zeros(n, np.uint64)
    for i in range(K):
        x = bb[i:i + n] & np.uint64(3)
        f = (f << np.uint64(2)) | x
        r = r | ((np.uint64(3) - x) << np.uint64(2 * i))
    c = np.cumsum(np.concatenate([[0], bad.astype(np.int64)]))
    v = (c[K:K + n] - c[:n]) == 0
    return f, r, v


def read_fasta(path):
    name, buf = None, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                if name is not None:
                    yield name, ''.join(buf)
                name, buf = line[1:].rstrip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, ''.join(buf)


def panel_a(ens, gene_id='ENSMUSG00000048856', R=101):
    import pysam
    W = R - K + 1
    gi = int(np.nonzero(GID == gene_id)[0][0]); G = GENES[gi]
    # TSL1 protein-coding transcripts and exons (e9rs.parse_gtf rule)
    pat = re.compile(r'(\w+) "([^"]*)"'); tx = {}
    with gzip.open(f'{ens}/Mus_musculus.GRCm38.102.gtf.gz', 'rt') as fh:
        for line in fh:
            if line[0] == '#' or gene_id not in line:
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] != 'exon':
                continue
            a = dict(pat.findall(f[8]))
            t = tx.setdefault(a['transcript_id'], dict(chrom=f[0], strand=f[6], exons=[],
                              keep=(a.get('gene_biotype') == 'protein_coding' and a.get('transcript_biotype') == 'protein_coding'
                                    and a.get('transcript_support_level', 'NA').split()[0] == '1')))
            t['exons'].append((int(f[3]), int(f[4])))
    seqs, gene_of = {}, {}
    for fn in ('Mus_musculus.GRCm38.cdna.all.fa.gz', 'Mus_musculus.GRCm38.ncrna.fa.gz'):
        for name, s in read_fasta(f'{ens}/{fn}'):
            tid = name.split()[0].split('.')[0]; mm = re.search(r'gene:(\S+)', name)
            gene_of[tid] = mm.group(1).split('.')[0] if mm else tid; seqs[tid] = s
    keep = [t for t, v in tx.items() if v['keep'] and t in seqs and len(seqs[t]) == sum(e - s + 1 for s, e in v['exons'])]
    keep = sorted(keep, key=lambda t: -len(seqs[t]))
    assert [len(seqs[t]) for t in keep] == G['tx_lens'] and keep[0] == G['ref_tx'], (keep, G)
    ref = keep[0]; rs = seqs[ref]; L = len(rs); b = LUT[np.frombuffer(rs.encode(), np.uint8)]
    f, r, v = codes(b); can = np.minimum(f, r); n = len(f); pos = np.arange(n)
    constit = v.copy()
    for t in keep[1:]:
        f2, _, v2 = codes(LUT[np.frombuffer(seqs[t].encode(), np.uint8)])
        constit &= np.isin(f, f2[v2])
    margin = (pos >= 100) & (pos <= L - 100 - K)
    # transcriptome uniqueness: canonical k-mer of this gene present in any transcript of another gene.
    # Exact string matching (Aho-Corasick) on both strands is equivalent to the canonical 2-bit comparison of e9rs.py.
    import ahocorasick
    comp = str.maketrans('ACGT', 'TGCA')

    def automaton(idx):
        A = ahocorasick.Automaton()
        for i in idx:
            km = rs[i:i + K].upper(); A.add_word(km, int(can[i])); A.add_word(km.translate(comp)[::-1], int(can[i]))
        A.make_automaton(); return A
    A = automaton(np.nonzero(v)[0]); seen_other = set()
    for tid, sq in seqs.items():
        if gene_of[tid] != gene_id:
            for _, c in A.iter(sq.upper()):
                seen_other.add(c)
    shared = np.array([bool(v[i]) and int(can[i]) in seen_other for i in range(n)])
    cand = constit & margin & ~shared
    # genome copies (canonical) in the GRCm38 primary assembly, for candidate k-mers
    A = automaton(np.nonzero(cand)[0]); gcount = {}
    with gzip.open(f'{ens}/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz', 'rb') as fh:
        data = fh.read()
    for rec in data.split(b'>')[1:]:
        nl = rec.index(b'\n'); sq = rec[nl + 1:].replace(b'\n', b'').decode('ascii').upper()
        for _, c in A.iter(sq):
            gcount[c] = gcount.get(c, 0) + 1
    del data
    gc_pos = np.array([gcount.get(int(can[i]), 0) if cand[i] else 0 for i in range(n)], np.int64)
    genome_ok = cand & (gc_pos <= 1)
    # founder variants on the transcript's exons (e9rs.stream_variants rule: genotype not 0/0 or missing;
    # interval pos .. pos + len(REF) - 1, plus 1 for indels)
    t = tx[ref]; ex = sorted(t['exons'])
    gcoord = np.concatenate([np.arange(s, e + 1) for s, e in ex]); gcoord = gcoord if t['strand'] == '+' else gcoord[::-1]
    chrom = t['chrom']; lo, hi = int(gcoord.min()), int(gcoord.max())
    cov = np.zeros((7, L), bool); varrows = []
    for bi, (letter, strain) in enumerate(FOUNDERS):
        for kind, extra in (('snps.dbSNP142.vcf.gz', 0), ('indels.dbSNP142.normed.vcf.gz', 1)):
            vf = pysam.VariantFile(f'{MGP}/{strain}.mgp.v5.{kind}')
            for rec in vf.fetch(chrom, lo - 200, hi + 200):
                gt = list(rec.samples.values())[0]['GT']
                if gt is None or all(a is None for a in gt) or all(a == 0 for a in gt if a is not None):
                    continue
                a0, a1 = rec.pos, rec.pos + len(rec.ref) - 1 + extra
                hitp = np.nonzero((gcoord >= a0) & (gcoord <= a1))[0]
                if len(hitp):
                    cov[bi, hitp] = True
                    varrows.append((letter, strain, 'indel' if extra else 'SNP', rec.pos, int(hitp.min()), int(hitp.max()),
                                    rec.ref, ','.join(rec.alts or ()), '/'.join('.' if x is None else str(x) for x in gt),
                                    ','.join(rec.filter.keys())))
    mask = np.zeros(n, np.uint8)
    for bi in range(7):
        c = np.cumsum(np.concatenate([[0], cov[bi].astype(np.int64)]))
        mask |= (((c[K:K + n] - c[:n]) > 0).astype(np.uint8) << bi)

    def spaced(ok):
        keepm = np.zeros(n, bool); last = -10 ** 9
        for j in np.nonzero(ok)[0]:
            if j - last >= W:
                keepm[j] = True; last = j
        return keepm
    in_ref = spaced(genome_ok); in_vf = spaced(genome_ok & (mask == 0))
    # check against the stored design file
    _, P = gene_probe_counts(f'probes_s{W}.npz')
    s = P['gene'] == gi
    fz = pd.DataFrame({'pos': P['pos'][s], 'in_ref': P['in_ref'][s], 'in_vf': P['in_vf'][s], 'mask': P['mask'][s], 'kmer': P['kmer'][s]})
    mine_ref = np.nonzero(in_ref)[0]; mine_vf = np.nonzero(in_vf)[0]
    assert list(mine_ref) == list(fz.pos[fz.in_ref]), ('ref', list(mine_ref), list(fz.pos[fz.in_ref]))
    assert list(mine_vf) == list(fz.pos[fz.in_vf]), ('vf', list(mine_vf), list(fz.pos[fz.in_vf]))
    assert list(mask[mine_ref]) == list(fz['mask'][fz.in_ref]), 'mask mismatch'
    assert all(rs[p:p + K] == k for p, k in zip(fz.pos, fz.kmer)), 'kmer mismatch'
    print('a: rebuilt probes equal the stored design file:', len(mine_ref), 'ref-unique,', len(mine_vf), 'PLUCK probes')
    status = np.where(~margin, 'end margin', np.where(~constit, 'not in every isoform', np.where(shared, 'shared with another gene',
                      np.where(gc_pos > 1, 'more than one genome copy', np.where(mask > 0, 'overlaps a founder variant', 'eligible')))))
    pa = pd.DataFrame({'pos': pos, 'status': status, 'founder_mask': mask,
                       'founders': [''.join(l for i, (l, _) in enumerate(FOUNDERS) if mk >> i & 1) for mk in mask],
                       'genome_copies': gc_pos, 'ref_unique_probe': in_ref.astype(int), 'pluck_probe': in_vf.astype(int)})
    path = f'{OUT}/fig1a_gene_positions.tsv'
    header(path, f'rebuilt from Ensembl 102 + MGP v5 with the e9rs.py rules; checked equal to respacing/probes/probes_s{W}.npz',
           [f'gene {G["gene_name"]} {gene_id}, reference transcript {ref} (longest TSL1 isoform, {L} nt), design R = {R}, spacing {W}',
            'one row per 25-mer start on the reference transcript; status is the first rule the k-mer fails, in the order of e9rs.py',
            'genome_copies is filled only for k-mers that pass the first three rules; founders = founder letters whose variants the 25-mer overlaps',
            'founder letters: A A/J, C 129S1/SvImJ, D NOD/ShiLtJ, E NZO/HlLtJ, F CAST/EiJ, G PWK/PhJ, H WSB/EiJ (B = C57BL/6J = reference)'])
    pa.to_csv(path, sep='\t', index=False, mode='a')
    path = f'{OUT}/fig1a_gene_variants.tsv'
    header(path, f'{MGP}/<strain>.mgp.v5.{{snps.dbSNP142,indels.dbSNP142.normed}}.vcf.gz, tabix region {chrom}:{lo - 200}-{hi + 200}',
           ['tx_start/tx_end = 0-based positions on the reference transcript covered by the call (indels extended by 1, as in e9rs.py)',
            'every call with a genotype other than 0/0 or missing is kept, whatever its FILTER, as in e9rs.py'])
    pd.DataFrame(varrows, columns=['founder', 'strain', 'kind', 'genome_pos', 'tx_start', 'tx_end', 'ref', 'alt', 'gt', 'filter']
                 ).sort_values(['tx_start', 'founder']).to_csv(path, sep='\t', index=False, mode='a')
    # isoform structure on the reference transcript coordinate: which reference positions each isoform contains
    iso = []
    for tname in keep:
        te = sorted(tx[tname]['exons']); gco = np.concatenate([np.arange(s, e + 1) for s, e in te])
        inside = np.isin(gcoord, gco); edges = np.flatnonzero(np.diff(np.concatenate([[0], inside.astype(int), [0]])))
        for s0, e0 in zip(edges[::2], edges[1::2]):
            iso.append((tname, len(seqs[tname]), int(s0), int(e0)))
    ref_ex = []
    acc = 0
    for s0, e0 in (ex if t['strand'] == '+' else ex[::-1]):
        ref_ex.append((acc, acc + e0 - s0 + 1)); acc += e0 - s0 + 1
    path = f'{OUT}/fig1a_gene_isoforms.tsv'
    header(path, 'Ensembl 102 GTF exons of the TSL1 protein-coding isoforms', ['segment = [start, end) on the reference transcript '
           'that the isoform contains; exon rows give the reference transcript exon boundaries'])
    pd.DataFrame(iso + [('exon:' + ref, L, s0, e0) for s0, e0 in ref_ex], columns=['transcript', 'length', 'start', 'end']
                 ).to_csv(path, sep='\t', index=False, mode='a')
    print(pd.Series(status).value_counts().to_dict())


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--ens'); ap.add_argument('--skip-a', action='store_true')
    ap.add_argument('--only', choices=['a', 'c', 'd'])
    a = ap.parse_args()
    if a.only in (None, 'c'):
        panel_c()
    if a.only in (None, 'd'):
        panel_d()
    if a.only in (None, 'a') and not a.skip_a:
        panel_a(a.ens)
