#!/usr/bin/env python3
"""Probe design and per-sample estimation with gene probes that avoid founder variants, compared with reference-unique probes.

Data: Keele et al. 2020, GSE140873 / PRJNA591379 (CC lung, liver, kidney; 50 bp SE, reverse-stranded).
Truth: the published per-sample RSEM estimates on strain-specific pseudo-genomes (GEO *_gene_quant.txt.gz).

Two training-free estimators, identical except for the founder-variant filter:
  (a) "ref":  25-mers in every TSL1 protein-coding isoform of the gene, unique to the gene in the
              Ensembl 102 transcriptome (cDNA + ncRNA, both strands) and occurring at most once in the
              GRCm38 genome (Kmerator-style; no variant filter);
  (b) "vf":   the subset of (a) overlapping no SNP or indel of any CC founder (Sanger MGP v5).
Each set is thinned to spaced probes (>= SPACING nt apart, >= 100 nt from transcript ends, along the
longest TSL1 isoform). Per gene: est = l~/(m (R-k+1)) * sum of probe counts, l~ = mean(L - F + 1).
Probes are counted antisense to the mRNA (read orientation of a reverse-stranded library) in a
per-sample msBWT built with pymsbwt.

Scoring: per sample, genes are stratified by the strain's founder haplotype at the gene, inferred from
MegaMUGA/GigaMUGA genotypes of Srivastava et al. 2017 (qtl2data CC): "B6" (B6 or locally identical
to it), "nonB6", and "nonB6_hit" (nonB6 and at least one ref probe carries that founder's variant).

Usage (Python >= 3.14 for pymsbwt; numpy, scipy):
  python e9.py run --work DIR --srr SRR10523526 SRR10523478      # prep (once) + samples + score
  python e9.py run --work DIR --titles all --tissue liver          # every liver RNA sample
  python e9.py prep|sample|score ...                               # single steps
Downloads stay under --work; counts/results are small. FASTQ and msBWT are deleted after counting
unless --keep-bwt.
"""
import argparse, gzip, json, os, re, shutil, subprocess, sys, time, zipfile
import numpy as np

K = 25
SPACING = 101
SPACINGS = [26, 51, 77, 101, 126]
END_MARGIN = 100
CHROMS = [str(i) for i in range(1, 20)] + ['X']
ENS = 'https://ftp.ensembl.org/pub/release-102'
URLS = {
    'cdna.fa.gz': f'{ENS}/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz',
    'ncrna.fa.gz': f'{ENS}/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz',
    'genes.gtf.gz': f'{ENS}/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz',
    'dna.fa.gz': f'{ENS}/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz',
    'cc.zip': 'https://raw.githubusercontent.com/rqtl/qtl2data/main/CC/cc.zip',
    'geo_filelist.txt': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140873/suppl/filelist.txt',
    'runinfo.tsv': ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA591379'
                    '&result=read_run&fields=run_accession,sample_title,read_count,fastq_ftp&format=tsv'),
}
MGP = 'https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs'
# qtl2 founder letters; B = C57BL/6J = reference. Bit order of the variant mask follows this list.
FOUNDERS = [('A', 'A_J'), ('C', '129S1_SvImJ'), ('D', 'NOD_ShiLtJ'), ('E', 'NZO_HlLtJ'),
            ('F', 'CAST_EiJ'), ('G', 'PWK_PhJ'), ('H', 'WSB_EiJ')]
FBIT = {l: 1 << i for i, (l, _) in enumerate(FOUNDERS)}

LUT = np.full(256, 4, np.uint8)
for i, ch in enumerate(b'ACGT'):
    LUT[ch] = i
    LUT[ch + 32] = i
RC = bytes.maketrans(b'ACGTacgtN', b'TGCAtgcaN')


def log(w, **kw):
    kw['t'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    print(json.dumps(kw), flush=True)
    with open(os.path.join(w, 'log.jsonl'), 'a') as f:
        f.write(json.dumps(kw) + '\n')


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')


def fetch(url, path):
    if os.path.exists(path):
        return path
    for attempt in range(6):
        try:
            sh(f"curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors -C - -o '{path}.part' '{url}'")
            if path.endswith('.gz'):
                sh(f"gzip -t '{path}.part'")
            os.replace(path + '.part', path)
            return path
        except subprocess.CalledProcessError:
            if path.endswith('.gz') and os.path.exists(path + '.part') and attempt >= 2:
                os.remove(path + '.part')  # corrupt partial: start over
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f'download failed: {url}')


def seqcodes(b):
    """b: uint8 0..4 (4 = N). Returns forward code, reverse-complement code, valid mask per k-mer start."""
    n = len(b) - K + 1
    if n <= 0:
        e = np.empty(0, np.uint64)
        return e, e, np.empty(0, bool)
    f = np.zeros(n, np.uint64)
    r = np.zeros(n, np.uint64)
    bb = (b & 3).astype(np.uint64)
    for i in range(K):
        x = bb[i:i + n]
        f = (f << np.uint64(2)) | x
        r |= (np.uint64(3) - x) << np.uint64(2 * i)
    cs = np.concatenate([[0], np.cumsum(b == 4, dtype=np.int64)])
    return f, r, (cs[K:K + n] - cs[:n]) == 0


def enc(s):
    return LUT[np.frombuffer(s.encode() if isinstance(s, str) else s, np.uint8)]


def read_fasta(path):
    name, parts = None, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                if name:
                    yield name, ''.join(parts)
                name, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if name:
        yield name, ''.join(parts)


# ---------------------------------------------------------------- prep
def parse_gtf(path):
    tx = {}
    pat = re.compile(r'(\w+) "([^"]*)"')
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] != 'exon' or f[0] not in CHROMS:
                continue
            a = dict(pat.findall(f[8]))
            t = tx.get(a['transcript_id'])
            if t is None:
                tsl = a.get('transcript_support_level', 'NA')
                t = tx[a['transcript_id']] = dict(
                    gene_id=a['gene_id'], gene_name=a.get('gene_name', a['gene_id']), chrom=f[0], strand=f[6],
                    keep=(a.get('gene_biotype') == 'protein_coding' and a.get('transcript_biotype') == 'protein_coding'
                          and tsl.split()[0] == '1'), exons=[])
            t['exons'].append((int(f[3]), int(f[4])))
    return tx


def tx_genome_coords(t):
    ex = sorted(t['exons'])
    g = np.concatenate([np.arange(s, e + 1, dtype=np.int64) for s, e in ex])
    return g if t['strand'] == '+' else g[::-1]


def stream_variants(w, letter, strain):
    out = os.path.join(w, 'ref', f'var_{letter}.npy')
    if os.path.exists(out):
        return np.load(out)
    parts = []
    for kind, extra in (('snps.dbSNP142.vcf.gz', 0), ('indels.dbSNP142.normed.vcf.gz', 1)):
        url = f'{MGP}/{strain}.mgp.v5.{kind}'
        tmp = os.path.join(w, 'ref', f'var_{letter}_{extra}.txt')
        local = os.path.join(w, 'ref', 'vcf', f'{strain}.mgp.v5.{kind}')
        src = f"gzip -dc '{local}'" if os.path.exists(local) else f"curl -fsSL --retry 5 --retry-all-errors '{url}' | gzip -dc"
        awk = ('!/^#/ { split($10,g,":"); if (g[1]=="0/0" || g[1]=="./." || g[1]==".") next; c=$1;'
               ' if (c=="X") c=20; else if (c !~ /^[0-9]+$/) next; print c, $2, $2+length($4)-1+' + str(extra) + ' }')
        for attempt in range(4):
            try:
                sh(f"set -o pipefail; {src} | awk -F'\\t' '{awk}' > '{tmp}'")
                break
            except subprocess.CalledProcessError:
                if attempt == 3:
                    raise
                time.sleep(10)
        parts.append(np.fromfile(tmp, dtype=np.int64, sep=' ').reshape(-1, 3))
        os.remove(tmp)
    v = np.concatenate(parts).astype(np.int32)
    np.save(out, v)
    return v


def prep(w, a):
    ref = os.path.join(w, 'ref')
    os.makedirs(ref, exist_ok=True)
    if os.path.exists(os.path.join(ref, 'probes.npz')):
        return
    t0 = time.time()
    for fn, url in URLS.items():
        fetch(url, os.path.join(ref, fn))
    log(w, step='prep_download', s=round(time.time() - t0))

    # transcripts and genes
    tx = parse_gtf(os.path.join(ref, 'genes.gtf.gz'))
    seqs = {}
    gene_of = {}
    for fn in ('cdna.fa.gz', 'ncrna.fa.gz'):
        for name, s in read_fasta(os.path.join(ref, fn)):
            tid = name.split()[0].split('.')[0]
            m = re.search(r'gene:(\S+)', name)
            gene_of[tid] = m.group(1).split('.')[0] if m else tid
            seqs[tid] = s
    genes = {}
    for tid, t in tx.items():
        if t['keep'] and tid in seqs and len(seqs[tid]) == sum(e - s + 1 for s, e in t['exons']):
            genes.setdefault(t['gene_id'], []).append(tid)
    log(w, step='prep_annot', tsl1_pc_genes=len(genes))

    # transcriptome uniqueness: canonical 25-mers present in more than one gene (any biotype)
    per_gene = {}
    for tid, s in seqs.items():
        per_gene.setdefault(gene_of[tid], []).append(tid)
    allc = []
    for gid, tids in per_gene.items():
        cs = []
        for tid in tids:
            f, r, v = seqcodes(enc(seqs[tid]))
            cs.append(np.minimum(f, r)[v])
        allc.append(np.unique(np.concatenate(cs)))
    allc = np.sort(np.concatenate(allc))
    multi = np.unique(allc[1:][allc[1:] == allc[:-1]])
    del allc
    log(w, step='prep_tx_unique', multigene_kmers=int(len(multi)))

    # candidates per gene: constitutive (every TSL1 pc isoform), transcriptome-unique, away from ends
    G = []  # gene records
    cand_gene, cand_pos, cand_can = [], [], []
    for gid, tids in genes.items():
        tids = sorted(tids, key=lambda t: -len(seqs[t]))
        rs = seqs[tids[0]]
        b = enc(rs)
        f, r, v = seqcodes(b)
        ok = v.copy()
        for t in tids[1:]:
            f2, _, v2 = seqcodes(enc(seqs[t]))
            ok &= np.isin(f, f2[v2])
        L = len(rs)
        pos = np.arange(len(f))
        ok &= (pos >= END_MARGIN) & (pos <= L - END_MARGIN - K)
        can = np.minimum(f, r)
        j = np.minimum(np.searchsorted(multi, can), len(multi) - 1)
        ok &= multi[j] != can  # multi is sorted; np.isin would re-sort it per gene
        idx = np.nonzero(ok)[0]
        t = tx[tids[0]]
        gc = tx_genome_coords(t)
        G.append(dict(gene_id=gid, gene_name=t['gene_name'], chrom=t['chrom'], strand=t['strand'],
                      mid=int((gc.min() + gc.max()) // 2), ref_tx=tids[0], tx_lens=[len(seqs[x]) for x in tids]))
        cand_gene.append(np.full(len(idx), len(G) - 1, np.int32))
        cand_pos.append(idx.astype(np.int32))
        cand_can.append(can[idx])
    cand_gene = np.concatenate(cand_gene); cand_pos = np.concatenate(cand_pos); cand_can = np.concatenate(cand_can)
    del multi
    log(w, step='prep_candidates', n=int(len(cand_can)))

    # genome uniqueness: occurrences of each candidate (canonical) in GRCm38 primary assembly
    U, inv = np.unique(cand_can, return_inverse=True)
    gcount = np.zeros(len(U), np.int64)
    chrom_len = {}
    dna = os.path.join(ref, 'dna.fa')
    if not os.path.exists(dna):
        sh(f"gzip -dc '{dna}.gz' > '{dna}.part' && mv '{dna}.part' '{dna}'")
    with open(dna, 'rb') as fh:
        data = fh.read()
    starts = [m.start() for m in re.finditer(rb'>', data)] + [len(data)]
    for i in range(len(starts) - 1):
        rec = data[starts[i]:starts[i + 1]]
        nl = rec.index(b'\n')
        name = rec[1:nl].split()[0].decode()
        b = LUT[np.frombuffer(rec[nl + 1:].replace(b'\n', b''), np.uint8)]
        chrom_len[name] = len(b)
        step = 20_000_000
        for s in range(0, len(b), step):
            f, r, v = seqcodes(b[s:s + step + K - 1])
            c = np.minimum(f, r)[v]
            j = np.searchsorted(U, c)
            j[j == len(U)] = 0
            hit = U[j] == c
            gcount += np.bincount(j[hit], minlength=len(U))
    del data
    genome_ok = gcount[inv] <= 1
    log(w, step='prep_genome_unique', pass_frac=round(float(genome_ok.mean()), 4))

    # founder variant mask per candidate (bit per founder: k-mer overlaps a SNP/indel of that founder)
    var = {l: stream_variants(w, l, s) for l, s in FOUNDERS}
    log(w, step='prep_variants', n={l: int(len(v)) for l, v in var.items()})
    mask = np.zeros(len(cand_can), np.uint8)
    gene_chrom = np.array([g['chrom'] for g in G])
    gb = np.searchsorted(cand_gene, np.arange(len(G) + 1))  # candidates are grouped by gene, in order
    for ci, ch in enumerate(CHROMS):
        cnum = 20 if ch == 'X' else int(ch)
        arr = np.zeros(chrom_len[ch] + 2, np.uint8)
        for l, v in var.items():
            vv = v[v[:, 0] == cnum]
            d = np.bincount(vv[:, 1], minlength=len(arr) + 1)[:len(arr) + 1].astype(np.int32)
            d -= np.bincount(np.minimum(vv[:, 2] + 1, len(arr)), minlength=len(arr) + 1)[:len(arr) + 1].astype(np.int32)
            arr |= (np.cumsum(d[:len(arr)]) > 0).astype(np.uint8) * np.uint8(FBIT[l])
        for gi in np.nonzero(gene_chrom == ch)[0]:
            sel = np.arange(gb[gi], gb[gi + 1])
            if not len(sel):
                continue
            m = arr[tx_genome_coords(tx[G[gi]['ref_tx']])]
            win = np.zeros(len(m) - K + 1, np.uint8)
            for k in range(K):
                win |= m[k:k + len(win)]
            mask[sel] = win[cand_pos[sel]]
        del arr

    # spaced selection for each set
    def spaced(okmask):
        keep = np.zeros(len(cand_can), bool)
        for gi in range(len(G)):
            sel = gb[gi] + np.nonzero(okmask[gb[gi]:gb[gi + 1]])[0]
            last = -10 ** 9
            for j in sel:
                if cand_pos[j] - last >= SPACING:
                    keep[j] = True
                    last = cand_pos[j]
        return keep
    order = np.lexsort((cand_pos, cand_gene))
    cand_gene, cand_pos, cand_can, genome_ok, mask = (x[order] for x in (cand_gene, cand_pos, cand_can, genome_ok, mask))
    np.savez(os.path.join(ref, 'pool.npz'), gene=cand_gene, pos=cand_pos, can=cand_can, genome_ok=genome_ok, mask=mask)
    global SPACING
    for SPACING in SPACINGS:
        in_ref = spaced(genome_ok)
        in_vf = spaced(genome_ok & (mask == 0))
        keep = in_ref | in_vf
        kmers = np.array([seqs[G[g]['ref_tx']][p:p + K] for g, p in zip(cand_gene[keep], cand_pos[keep])])
        np.savez_compressed(os.path.join(ref, f'probes_s{SPACING}.npz'), gene=cand_gene[keep], pos=cand_pos[keep], mask=mask[keep],
                            in_ref=in_ref[keep], in_vf=in_vf[keep], kmer=kmers)
        log(w, step='prep_spacing', spacing=SPACING, ref_probes=int(in_ref.sum()), vf_probes=int(in_vf.sum()),
            genes_vf=int(len(np.unique(cand_gene[in_vf]))))
    open(os.path.join(ref, 'probes.npz'), 'w').close()
    with open(os.path.join(ref, 'genes.json'), 'w') as f:
        json.dump(G, f)
    log(w, step='prep_done', s=round(time.time() - t0), genes=len(G),
        ref_probes=int(in_ref.sum()), vf_probes=int(in_vf.sum()),
        )
    os.remove(dna)


def fetch_all(w, srrs):
    ref = os.path.join(w, 'ref')
    for fn, url in URLS.items():
        fetch(url, os.path.join(ref, fn))
    os.makedirs(os.path.join(ref, 'vcf'), exist_ok=True)
    for _, strain in FOUNDERS:
        for kind in ('snps.dbSNP142.vcf.gz', 'indels.dbSNP142.normed.vcf.gz'):
            fetch(f'{MGP}/{strain}.mgp.v5.{kind}', os.path.join(ref, 'vcf', f'{strain}.mgp.v5.{kind}'))
    os.makedirs(os.path.join(w, 'fastq'), exist_ok=True)
    info = runinfo(w)
    for s in srrs:
        fetch('https://' + info[s]['fastq_ftp'].split(';')[0], os.path.join(w, 'fastq', f'{s}.fastq.gz'))
        strain, tissue, _ = info[s]['sample_title'].split('_')
        truth_path(w, strain, tissue)
    log(w, step='fetch_done', samples=len(srrs))


def truth_path(w, strain, tissue):
    geo = {}
    for line in open(os.path.join(w, 'ref', 'geo_filelist.txt')):
        m = re.search(r'(GSM\d+)_(CC\d+)_(\w+?)_gene_quant\.txt\.gz', line)
        if m:
            geo[(m.group(2), m.group(3))] = m.group(0)
    gfile = geo[(strain, tissue)]
    os.makedirs(os.path.join(w, 'truth'), exist_ok=True)
    return fetch(f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gfile[:7]}nnn/{gfile.split('_')[0]}/suppl/{gfile}",
                 os.path.join(w, 'truth', gfile))


# ---------------------------------------------------------------- sample
def runinfo(w):
    rows = {}
    with open(os.path.join(w, 'ref', 'runinfo.tsv')) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            r = dict(zip(hdr, line.rstrip('\n').split('\t')))
            if r['sample_title'].endswith('_rna'):
                rows[r['run_accession']] = r
    return rows


def sample(w, srr, a):
    out = os.path.join(w, 'counts', f'{srr}.npz')
    if os.path.exists(out):
        return
    os.makedirs(os.path.join(w, 'counts'), exist_ok=True)
    tmp = os.path.join(w, 'tmp', srr)
    os.makedirs(tmp, exist_ok=True)
    info = runinfo(w)[srr]
    bwt = os.path.join(tmp, 'bwt')
    t0 = t1 = time.time()
    rlen = None
    if not os.path.exists(os.path.join(bwt, 'comp_msbwt.npy')):  # reuse a kept msBWT
        pre = os.path.join(w, 'fastq', f'{srr}.fastq.gz')
        gz = pre if os.path.exists(pre) else fetch('https://' + info['fastq_ftp'].split(';')[0], os.path.join(tmp, 'r.fastq.gz'))
        t1 = time.time()
        fq = os.path.join(tmp, 'r.fq')
        sh(f"gzip -dc '{gz}' > '{fq}'")
        with open(fq) as f:
            f.readline()
            rlen = len(f.readline().strip())
        if gz != pre:
            os.remove(gz)
        exe = shutil.which('msbwt') or os.path.join(os.path.dirname(sys.executable), 'msbwt')
        sh(f"'{exe}' cffq -p {a.threads} -u -c '{bwt}' '{fq}' > '{tmp}/build.log' 2>&1")
        os.remove(fq)
    t2 = time.time()
    import MUSCython.MultiStringBWTCython as MSBWT
    msbwt = MSBWT.loadBWT(bwt)
    P = dict(np.load(os.path.join(w, 'ref', 'probes.npz')))  # npz re-reads on every key access
    q = [k.translate(RC)[::-1] if a.orientation == 'antisense' else k for k in P['kmer'].tolist()]
    counts = np.array([msbwt.countOccurrencesOfSeq(s) for s in q], np.int64)
    t3 = time.time()
    if rlen is None:  # uniform reads: symbols per read minus the '$' terminator
        tc = np.load(os.path.join(bwt, 'totalCounts.npy'))
        rlen = int(round(tc.sum() / tc[0])) - 1
    np.savez_compressed(out, counts=counts, read_len=rlen, title=info['sample_title'], reads=int(info['read_count']))
    log(w, step='sample', srr=srr, title=info['sample_title'], read_len=rlen, download_s=round(t1 - t0),
        build_s=round(t2 - t1), query_s=round(t3 - t2), n_queries=len(q))
    if not a.keep_bwt:
        shutil.rmtree(tmp)


# ---------------------------------------------------------------- score
def load_cc(w):
    z = zipfile.ZipFile(os.path.join(w, 'ref', 'cc.zip'))
    rd = lambda n: [l.split(',') for l in z.read(n).decode().splitlines() if not l.startswith('#')]
    cc = {}
    for ch in CHROMS:
        g = rd(f'cc_geno{ch}.csv'); fg = rd(f'MMnGM/MMnGM_foundergeno{ch}.csv'); pm = rd(f'MMnGM/MMnGM_pmap{ch}.csv')
        fdict = {r[0]: r[1:] for r in fg[1:]}
        pdict = {r[0]: float(r[2]) for r in pm[1:]}
        rows = [r for r in g[1:] if r[0] in fdict and r[0] in pdict and len(set(fdict[r[0]]) - {'-'}) > 1]
        rows.sort(key=lambda r: pdict[r[0]])
        cc[ch] = dict(strains=[s.split('/')[0] for s in g[0][1:]], pos=np.array([pdict[r[0]] for r in rows]),
                      geno=np.array([r[1:] for r in rows]), fgeno=np.array([fdict[r[0]] for r in rows]))
    return cc


def haplotype(cc, strain, ch, mid_bp, nmark=10):
    """Return (stratum, best founder letters) for a CC strain at a position."""
    d = cc[ch]
    if strain not in d['strains']:
        return 'unknown', ''
    key = ('_ok', strain)
    if key not in d:  # cache the strain's genotyped markers on this chromosome
        col = d['geno'][:, d['strains'].index(strain)]
        d[key] = (col, np.nonzero(col != '-')[0])
    col, ok = d[key]
    if len(ok) < nmark:
        return 'ambiguous', ''
    i = np.searchsorted(d['pos'][ok], mid_bp / 1e6)
    win = ok[max(0, i - nmark):i + nmark]
    near = win[np.argsort(np.abs(d['pos'][win] - mid_bp / 1e6))[:nmark]]
    fg = d['fgeno'][near]
    score = ((fg == col[near, None]) + 0.5 * (fg == '-')).sum(0)
    best = ''.join('ABCDEFGH'[i] for i in np.nonzero(score == score.max())[0])
    if score.max() < nmark - 2:
        return 'ambiguous', best
    if 'B' in best:
        return 'B6', best
    return ('nonB6', best) if score[1] <= score.max() - 3 else ('ambiguous', best)


def spearman(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(x, y).statistic) if len(x) > 2 else float('nan')


def score(w, a):
    P = dict(np.load(os.path.join(w, 'ref', 'probes.npz')))  # npz re-reads on every key access
    G = json.load(open(os.path.join(w, 'ref', 'genes.json')))
    cc = load_cc(w)
    bnd = np.searchsorted(P['gene'], np.arange(len(G) + 1))  # probes are sorted by gene
    names = np.array([g['gene_name'] for g in G])
    dup = {n for n, c in zip(*np.unique(names, return_counts=True)) if c > 1}
    rows, pergene = [], []
    for fn in sorted(os.listdir(os.path.join(w, 'counts'))):
        Z = np.load(os.path.join(w, 'counts', fn))
        srr = fn[:-4]
        strain, tissue, _ = str(Z['title']).split('_')
        R = int(Z['read_len']); F = a.frag_len or R
        tp = truth_path(w, strain, tissue)
        truth = {}
        with gzip.open(tp, 'rt') as f:
            f.readline()
            for line in f:
                x = line.split('\t')
                truth[x[0]] = float(x[4])
        c = Z['counts'].astype(float)
        est = {}
        for s in ('ref', 'vf'):
            sel = P['in_' + s]
            S = np.bincount(P['gene'][sel], weights=c[sel], minlength=len(G))
            m = np.bincount(P['gene'][sel], minlength=len(G))
            lt = np.array([np.mean([max(L - F + 1, K) for L in g['tx_lens']]) for g in G])
            with np.errstate(invalid='ignore', divide='ignore'):
                est[s] = np.where(m > 0, lt / (np.maximum(m, 1) * (R - K + 1)) * S, np.nan)
        recs = []
        for gi, g in enumerate(G):
            n = g['gene_name']
            if n in dup or n not in truth or np.isnan(est['ref'][gi]) or np.isnan(est['vf'][gi]):
                continue
            st, best = haplotype(cc, strain, g['chrom'], g['mid'])
            # hit: some ref probe of the gene overlaps a variant carried by every inferred founder
            hb = 0
            for l in best:
                hb |= FBIT.get(l, 0)
            mk = P['mask'][bnd[gi]:bnd[gi + 1]][P['in_ref'][bnd[gi]:bnd[gi + 1]]]
            hit = st == 'nonB6' and hb > 0 and bool(np.any((mk & hb) == hb))
            recs.append((n, st, best, hit, truth[n], est['ref'][gi], est['vf'][gi]))
        tr = np.array([r[4] for r in recs]); er = np.array([r[5] for r in recs]); ev = np.array([r[6] for r in recs])
        stv = np.array([r[1] for r in recs]); hitv = np.array([r[3] for r in recs])
        expr = tr >= a.min_truth
        for s, e in (('ref', er), ('vf', ev)):
            base = expr & (stv == 'B6') & (e > 0)
            C = float(np.median(tr[base] / e[base]))
            lr = np.log2((C * e + 1) / (tr + 1))
            for label, sm in (('B6', stv == 'B6'), ('nonB6', stv == 'nonB6'),
                              ('nonB6_hit', hitv), ('nonB6_nohit', (stv == 'nonB6') & ~hitv), ('all', stv != 'x')):
                q = expr & sm
                rows.append(dict(srr=srr, sample=str(Z['title']), set=s, stratum=label, n=int(q.sum()), C=round(C, 3),
                                 spearman=round(spearman(e[q], tr[q]), 4),
                                 med_log2=round(float(np.median(lr[q])), 4) if q.any() else None,
                                 frac_off2x=round(float(np.mean(np.abs(lr[q]) > 1)), 4) if q.any() else None))
        for r, e1, e2 in zip(recs, er, ev):
            pergene.append((srr, tissue, strain) + r[:4] + (r[4], e1, e2))
    os.makedirs(os.path.join(w, 'results'), exist_ok=True)
    with open(os.path.join(w, 'results', 'summary.tsv'), 'w') as f:
        keys = list(rows[0])
        f.write('\t'.join(keys) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[k]) for k in keys) + '\n')
    with open(os.path.join(w, 'results', 'per_gene.tsv'), 'w') as f:
        f.write('srr\ttissue\tstrain\tgene\tstratum\tfounders\thit\ttruth\test_ref\test_vf\n')
        for p in pergene:
            f.write('\t'.join(str(x) for x in p) + '\n')
    cohort(w, pergene, a)
    for r in rows:
        print('\t'.join(str(v) for v in r.values()))


def cohort(w, pergene, a):
    """Within-gene, across-strain accuracy (needs many samples of one tissue)."""
    out = []
    by = {}
    for p in pergene:
        by.setdefault((p[1], p[3]), []).append(p)
    for tissue in sorted({k[0] for k in by}):
        res = {('ref', 'mixed'): [], ('vf', 'mixed'): [], ('ref', 'allB6'): [], ('vf', 'allB6'): []}
        for (t, gname), ps in by.items():
            if t != tissue or len(ps) < a.min_samples:
                continue
            tr = np.log2(np.array([p[7] for p in ps]) + 1)
            if np.median([p[7] for p in ps]) < a.min_truth or tr.std() == 0:
                continue
            nb = sum(p[4] == 'nonB6' for p in ps)
            grp = 'mixed' if nb >= 3 and sum(p[4] == 'B6' for p in ps) >= 3 else ('allB6' if nb == 0 else None)
            if grp is None:
                continue
            for s, j in (('ref', 8), ('vf', 9)):
                e = np.log2(np.array([p[j] for p in ps]) + 1)
                if e.std() > 0:
                    res[(s, grp)].append(np.corrcoef(e, tr)[0, 1] ** 2)
        for (s, grp), v in res.items():
            if v:
                out.append(f'{tissue}\t{s}\t{grp}\t{len(v)}\t{np.median(v):.4f}')
    if out:
        with open(os.path.join(w, 'results', 'cohort.tsv'), 'w') as f:
            f.write('tissue\tset\tgene_group\tn_genes\tmedian_within_gene_r2\n' + '\n'.join(out) + '\n')
        print('\n'.join(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['fetch', 'list', 'prep', 'sample', 'score', 'run'])
    ap.add_argument('--work', required=True)
    ap.add_argument('--srr', nargs='*', default=[])
    ap.add_argument('--titles', help="'all' = every RNA run in PRJNA591379 (optionally with --tissue)")
    ap.add_argument('--tissue')
    ap.add_argument('--threads', type=int, default=os.cpu_count())
    ap.add_argument('--orientation', choices=['antisense', 'sense'], default='antisense')
    ap.add_argument('--frag-len', type=int, default=0, help='for effective length; default = read length (SE)')
    ap.add_argument('--min-truth', type=float, default=50)
    ap.add_argument('--min-samples', type=int, default=10)
    ap.add_argument('--keep-bwt', action='store_true')
    ap.add_argument('--task-index', type=int, help='run only the Nth (0-based) selected sample, e.g. $SLURM_ARRAY_TASK_ID')
    a = ap.parse_args()
    w = os.path.abspath(a.work)
    os.makedirs(w, exist_ok=True)
    os.makedirs(os.path.join(w, 'ref'), exist_ok=True)
    fetch(URLS['runinfo.tsv'], os.path.join(w, 'ref', 'runinfo.tsv'))
    srrs = list(a.srr)
    if a.titles == 'all':
        srrs += [k for k, r in sorted(runinfo(w).items(), key=lambda x: x[1]['sample_title'])
                 if not a.tissue or f'_{a.tissue}_' in r['sample_title']]
    if a.task_index is not None:
        srrs = [srrs[a.task_index]]
    if a.cmd == 'list':
        print(len(srrs))
        return
    if a.cmd == 'fetch':  # everything that needs internet, for clusters whose compute nodes are offline
        fetch_all(w, srrs)
        return
    if a.cmd in ('prep', 'run'):
        prep(w, a)
    if a.cmd in ('sample', 'run'):
        for s in srrs:
            sample(w, s, a)
    if a.cmd in ('score', 'run'):
        score(w, a)


if __name__ == '__main__':
    main()
