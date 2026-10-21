#!/usr/bin/env python3
"""Founder-diagnostic k-mers -> per-gene haplotype calls and sample-identity QC.

For each exonic SNP that is homozygous-alt (PASS) in a known subset of the 7 non-B6 CC founders
(MGP v5), the 25-mer centred on the SNP is built for the reference and the alt allele.
Counting both k-mers in a sample gives per-gene founder likelihoods; the founders that fit
best are the local haplotype call. Calls are scored against marker haplotypes (qtl2data CC,
inferred with e9.haplotype), and the genome-wide calls are matched against every CC strain.

Counts are exact k-mer occurrence counts in the reads (both orientations), i.e. what an
msBWT query of the same k-mers returns; here a streaming numpy counter is used instead of
building the index.

Usage: dk.py prep W ; dk.py count W SRR... ; dk.py score W SRR...
W must contain ref/{genes.gtf.gz,cdna.fa.gz,ncrna.fa.gz,chr*.fa.gz,cc.zip,runinfo.tsv},
vcf/<strain>.tsv (chrom pos ref alt gt filter, from getvcf.py) and fastq/<SRR>.fastq.gz.
"""
import gzip, json, os, re, subprocess, sys, time
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'e9'))
import e9  # reuse load_cc / haplotype only (read-only)

K, HALF = 25, 12
CHR = [str(i) for i in range(1, 20)] + ['X']  # genome-wide: autosomes 1-19 and X
LET = 'ABCDEFGH'  # B = C57BL/6J = reference
STRAIN = {'A': 'A_J', 'C': '129S1_SvImJ', 'D': 'NOD_ShiLtJ', 'E': 'NZO_HlLtJ',
          'F': 'CAST_EiJ', 'G': 'PWK_PhJ', 'H': 'WSB_EiJ'}
LUT = np.full(256, 4, np.uint8)
for i, ch in enumerate(b'ACGT'):
    LUT[ch] = i
    LUT[ch + 32] = i
EPS, LLR = 0.02, np.log(100)


def kmers(c):
    """Canonical 2-bit k-mers along the last axis of uint8 codes (4 = N); returns (vals, valid)."""
    n = c.shape[-1] - K + 1
    f = np.zeros(c.shape[:-1] + (n,), np.uint64)
    r = np.zeros_like(f)
    for i in range(K):
        x = c[..., i:i + n].astype(np.uint64) & np.uint64(3)
        f = (f << np.uint64(2)) | x
        r |= (np.uint64(3) - x) << np.uint64(2 * i)
    bad = np.cumsum(c == 4, axis=-1, dtype=np.int32)
    bad = np.concatenate([np.zeros(c.shape[:-1] + (1,), np.int32), bad], -1)
    valid = (bad[..., K:] - bad[..., :n]) == 0
    return np.minimum(f, r), valid


def read_fasta(path):
    name, buf = None, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                if name:
                    yield name, ''.join(buf)
                name, buf = line[1:].rstrip(), []
            else:
                buf.append(line.rstrip())
    if name:
        yield name, ''.join(buf)


def log(w, **kw):
    kw['t'] = time.strftime('%H:%M:%S')
    print(json.dumps(kw), flush=True)
    with open(os.path.join(w, 'log.jsonl'), 'a') as f:
        f.write(json.dumps(kw) + '\n')


# ---------------------------------------------------------------- prep
def prep(w):
    genes, exons = {}, {}
    with gzip.open(os.path.join(w, 'ref', 'genes.gtf.gz'), 'rt') as fh:
        for line in fh:
            if line[0] == '#':
                continue
            x = line.split('\t')
            if x[2] != 'exon' or x[0] not in CHR or 'transcript_biotype "protein_coding"' not in x[8]:
                continue
            gid = re.search(r'gene_id "([^"]+)"', x[8]).group(1)
            name = re.search(r'gene_name "([^"]+)"', x[8]).group(1)
            genes[gid] = (name, x[0], x[6])
            exons.setdefault(gid, []).append((int(x[3]), int(x[4])))
    gl = sorted(genes)
    gidx = {g: i for i, g in enumerate(gl)}
    rows = []  # merged exon intervals: chrom, start, end, gene index
    for g in gl:
        iv = sorted(exons[g])
        cur = list(iv[0])
        for s, e in iv[1:]:
            if s <= cur[1] + 1:
                cur[1] = max(cur[1], e)
            else:
                rows.append((genes[g][1], cur[0], cur[1], gidx[g])); cur = [s, e]
        rows.append((genes[g][1], cur[0], cur[1], gidx[g]))
    mid = {g: (min(s for s, _ in exons[g]) + max(e for _, e in exons[g])) // 2 for g in gl}

    sites = []
    for ch in CHR:
        seq = np.frombuffer(next(read_fasta(os.path.join(w, 'ref', f'chr{ch}.fa.gz')))[1].upper().encode(), np.uint8)
        L = len(seq)
        parts = []  # pos, end, founder, snp_ok, alt
        for fi, l in enumerate(LET):
            if l == 'B':
                continue
            txt = os.path.join(w, 'vcf', f'{STRAIN[l]}.chr{ch}.txt')
            if not os.path.exists(txt):
                awk = ('$1=="%s" { ok = (length($3)==1 && length($4)==1 && $5=="1/1" && $6=="PASS"); '
                       'n = length($3); if (length($4) > n) n = length($4); '
                       'a = ok ? index("ACGT", $4) : 0; print $2, $2+n-1, %d, ok, a }' % (ch, fi))
                subprocess.run(f"awk -F'\\t' '{awk}' '{w}/vcf/{STRAIN[l]}.tsv' > '{txt}'", shell=True, check=True)
            parts.append(np.fromfile(txt, dtype=np.int64, sep=' ').reshape(-1, 5))
        R = np.concatenate(parts)
        R[:, 4] = np.where(R[:, 4] > 0, np.frombuffer(b'NACGT', np.uint8)[R[:, 4]], 0)
        allx = sorted((r[1], r[2]) for r in rows if r[0] == ch)
        ms, me = [], []
        for a_, b_ in allx:
            if ms and a_ <= me[-1] + 1:
                me[-1] = max(me[-1], b_)
            else:
                ms.append(a_); me.append(b_)
        ms, me = np.array(ms), np.array(me)
        jj = np.searchsorted(ms, R[:, 0], 'right') - 1
        exonic = (jj >= 0) & (R[:, 0] <= me[np.maximum(jj, 0)])
        # occupancy of any variant (any founder, any call) along the chromosome
        d = (np.bincount(np.clip(R[:, 0], 0, L), minlength=L + 2)
             - np.bincount(np.clip(R[:, 1] + 1, 0, L + 1), minlength=L + 2))
        occ = (np.cumsum(d) > 0).astype(np.int32)
        cs = np.concatenate([[0], np.cumsum(occ)])  # cs[p] = occupied positions before p
        # group records by position
        R = R[exonic][np.lexsort((R[exonic, 2], R[exonic, 0]))]
        up, start = np.unique(R[:, 0], return_index=True)
        stop = np.append(start[1:], len(R))
        ex = sorted(r for r in rows if r[0] == ch)
        exs = np.array([r[1] for r in ex]); exe = np.array([r[2] for r in ex]); exg = np.array([r[3] for r in ex])
        for p, a, b in zip(up, start, stop):
            G = R[a:b]
            if not G[:, 3].all() or len(set(G[:, 4])) != 1 or (G[:, 1] != p).any():
                continue
            if p - HALF < 1 or p + HALF > L:
                continue
            if cs[p + HALF + 1] - cs[p - HALF] != 1:  # another variant inside the k-mer window
                continue
            j = np.searchsorted(exs, p, 'right') - 1
            hits = [jj for jj in range(max(0, j - 50), j + 1) if exs[jj] <= p - HALF and exe[jj] >= p + HALF]
            allhit = {exg[jj] for jj in range(max(0, j - 50), j + 1) if exs[jj] <= p <= exe[jj]}
            if not hits or len(allhit) != 1:
                continue
            refk = seq[p - 1 - HALF:p + HALF]
            if LUT[refk].max() > 3:
                continue
            altk = refk.copy(); altk[HALF] = G[0, 4]
            sdp = int(sum(1 << int(f) for f in G[:, 2]))
            sites.append((exg[hits[0]], int(ch) if ch != 'X' else 20, p, sdp, refk.tobytes().decode(), altk.tobytes().decode()))
        log(w, step='sites', chrom=ch, n=len(sites))
    gi = np.array([s[0] for s in sites]); pos = np.array([s[2] for s in sites]); sdp = np.array([s[3] for s in sites])
    refc = kmers(LUT[np.frombuffer(''.join(s[4] for s in sites).encode(), np.uint8).reshape(-1, K)])[0][:, 0]
    altc = kmers(LUT[np.frombuffer(''.join(s[5] for s in sites).encode(), np.uint8).reshape(-1, K)])[0][:, 0]

    # transcriptome uniqueness: ref k-mer only in its own gene, alt k-mer nowhere
    q = np.unique(np.concatenate([refc, altc]))
    hitgene = {}
    tgene = {g: i for i, g in enumerate(gl)}
    for fa in ('cdna.fa.gz', 'ncrna.fa.gz'):
        buf, gids = [], []
        def flush():
            c = LUT[np.frombuffer(b''.join(buf), np.uint8)]
            v, ok = kmers(c)
            ga = np.repeat(np.array(gids), [len(b) for b in buf])[:len(v)]
            ii = np.searchsorted(q, v); ii[ii == len(q)] = 0
            m = ok & (q[ii] == v)
            for qi, g in zip(ii[m].tolist(), ga[m].tolist()):
                hitgene.setdefault(qi, set()).add(g)
        for name, s in read_fasta(os.path.join(w, 'ref', fa)):
            g = re.search(r'gene:(\S+?)(\.\d+)?\s', name + ' ').group(1)
            buf.append(s.encode() + b'N'); gids.append(tgene.get(g, -1 - len(gids)))
            if sum(len(b) for b in buf) > 2e7:
                flush(); buf, gids = [], []
        if buf:
            flush()
    ri = np.searchsorted(q, refc); ai = np.searchsorted(q, altc)
    keep = np.array([hitgene.get(r, set()) <= {g} and a not in hitgene and refc[k] != altc[k]
                     for k, (r, a, g) in enumerate(zip(ri.tolist(), ai.tolist(), gi.tolist()))])
    np.savez(os.path.join(w, 'sites.npz'), gene=gi[keep], pos=pos[keep], sdp=sdp[keep], ref=refc[keep], alt=altc[keep],
             chrom=np.array([s[1] for s in sites])[keep], genes=np.array(gl), names=np.array([genes[g][0] for g in gl]),
             gchrom=np.array([genes[g][1] for g in gl]), gmid=np.array([mid[g] for g in gl]))
    log(w, step='prep_done', candidates=len(sites), unique=int(keep.sum()), genes=int(len(set(gi[keep].tolist()))))


# ---------------------------------------------------------------- count
def count(w, srr):
    out = os.path.join(w, 'counts', f'{srr}.npz')
    if os.path.exists(out):
        return
    os.makedirs(os.path.dirname(out), exist_ok=True)
    S = np.load(os.path.join(w, 'sites.npz'))
    q = np.unique(np.concatenate([S['ref'], S['alt']]))
    tot = np.zeros(len(q), np.int64)
    t0 = time.time(); nreads = 0; L = None
    p = subprocess.Popen(f"gzip -dc '{w}/fastq/{srr}.fastq.gz' | awk 'NR%4==2'", shell=True, stdout=subprocess.PIPE, bufsize=1 << 24)
    while True:
        lines = p.stdout.readlines(1 << 25)
        if not lines:
            break
        L = L or max(len(x) for x in lines[:1000]) - 1
        c = LUT[np.frombuffer(b''.join(x[:-1].ljust(L, b'N')[:L] for x in lines), np.uint8).reshape(-1, L)]
        v, ok = kmers(c)
        v = v[ok]
        ii = np.searchsorted(q, v); ii[ii == len(q)] = 0
        tot += np.bincount(ii[q[ii] == v], minlength=len(q))
        nreads += len(lines)
    p.wait()
    np.savez(out, ref=tot[np.searchsorted(q, S['ref'])], alt=tot[np.searchsorted(q, S['alt'])], reads=nreads, read_len=L)
    log(w, step='count', srr=srr, reads=nreads, read_len=L, s=round(time.time() - t0))


# ---------------------------------------------------------------- score
def gene_ll(S, C):
    """Per-gene founder log-likelihoods (genes x 8) and per-gene diagnostic coverage."""
    carry = ((S['sdp'][:, None] >> np.arange(8)) & 1).astype(bool)  # site x founder carries alt
    pa = np.where(carry, 1 - EPS, EPS)
    ll_site = C['alt'][:, None] * np.log(pa) + C['ref'][:, None] * np.log(1 - pa)
    ng = len(S['genes'])
    ll = np.zeros((ng, 8)); np.add.at(ll, S['gene'], ll_site)
    cov = np.bincount(S['gene'], C['ref'] + C['alt'], minlength=ng)
    nsite = np.bincount(S['gene'], (C['ref'] + C['alt']) > 0, minlength=ng)
    return ll, cov, nsite


def score(w, srrs, mincov=10):
    S = dict(np.load(os.path.join(w, 'sites.npz')))
    cc = e9.load_cc(w)
    info = {}
    for line in open(os.path.join(w, 'ref', 'runinfo.tsv')):
        x = line.split('\t')
        info[x[0]] = x[1]
    strains = cc['1']['strains']
    rng = np.random.default_rng(1)
    res = {}
    for srr in srrs:
        C = np.load(os.path.join(w, 'counts', f'{srr}.npz'))
        strain = info[srr].split('_')[0]
        ll, cov, nsite = gene_ll(S, C)
        call = (ll >= ll.max(1, keepdims=True) - LLR)  # genes x founders within LLR of the best
        inf = (cov >= mincov) & (call.sum(1) < 8)
        rows, per = [], []
        for g in np.nonzero(inf)[0]:
            ch, mid = str(S['gchrom'][g]), int(S['gmid'][g])
            stratum, best = e9.haplotype(cc, strain, ch, mid)
            cs = ''.join(LET[i] for i in np.nonzero(call[g])[0])
            per.append((S['names'][g], ch, mid, int(cov[g]), int(nsite[g]), cs, stratum, best))
        tab = [r for r in per if r[6] in ('B6', 'nonB6')]
        consistent = np.array([bool(set(r[5]) & set(r[7])) for r in tab])
        single = np.array([len(r[5]) == 1 for r in tab])
        res[srr] = dict(strain=strain, reads=int(C['reads']), genes_called=len(per), genes_with_marker_truth=len(tab),
                        consistent=round(float(consistent.mean()), 3),
                        single_founder_calls=int(single.sum()),
                        single_founder_accuracy=round(float(consistent[single].mean()), 3),
                        median_call_set_size=float(np.median([len(r[5]) for r in tab])))
        for lo, hi in ((10, 30), (30, 100), (100, 10**9)):
            m = np.array([lo <= r[3] < hi for r in tab])
            res[srr][f'consistent_cov{lo}'] = [int(m.sum()), round(float(consistent[m].mean()), 3) if m.any() else None]
        # sample identity: fraction of called genes whose call fits each CC strain's marker haplotype
        trial = [r for r in per if len(r[5]) <= 4]
        fit = np.zeros((len(strains), len(trial)), bool); known = np.zeros_like(fit)
        for si, s in enumerate(strains):
            for gi_, r in enumerate(trial):
                _, best = e9.haplotype(cc, s, r[1], r[2])
                if best:
                    known[si, gi_] = True
                    fit[si, gi_] = bool(set(r[5]) & set(best))
        frac = fit.sum(1) / np.maximum(known.sum(1), 1)
        order = np.argsort(-frac)
        res[srr]['identity_top5'] = [(strains[i], round(float(frac[i]), 3)) for i in order[:5]]
        res[srr]['identity_rank_of_true'] = int(np.nonzero(np.array(strains)[order] == strain)[0][0]) + 1
        # genes needed: random subsets, how often the true strain is the unique top
        ti = strains.index(strain)
        need = {}
        for n in (10, 25, 50, 100, 200):
            if n > len(trial):
                break
            ok = 0
            for _ in range(200):
                sub = rng.choice(len(trial), n, replace=False)
                f = fit[:, sub].sum(1) / np.maximum(known[:, sub].sum(1), 1)
                ok += f[ti] > np.delete(f, ti).max()
            need[n] = ok / 200
        res[srr]['p_unique_top_by_ngenes'] = need
        with open(os.path.join(w, f'calls_{srr}.tsv'), 'w') as f:
            f.write('gene\tchrom\tmid\tcov\tsites_seen\tcall\tmarker_stratum\tmarker_best\n')
            for r in per:
                f.write('\t'.join(map(str, r)) + '\n')
    with open(os.path.join(w, 'summary.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    cmd, w = sys.argv[1], sys.argv[2]
    if cmd == 'prep':
        prep(w)
    elif cmd == 'count':
        for s in sys.argv[3:]:
            count(w, s)
    elif cmd == 'score':
        score(w, sys.argv[3:])
