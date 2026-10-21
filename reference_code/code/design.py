# Transcript-level sparse probe design (coordinate-hashed, vf) + union with 26 nt gene probes; compatibility map.
import sys, re, gzip, json, time, numpy as np
sys.path.insert(0, '.')
import e9rs
from e9rs import K, CHROMS, FOUNDERS, FBIT, seqcodes, enc, read_fasta, parse_gtf, tx_genome_coords
SP = 26
t0 = time.time()
tx = parse_gtf('ref/Mus_musculus.GRCm38.102.gtf.gz')
seqs, gene_of, btype = {}, {}, {}
for fn in ('ref/Mus_musculus.GRCm38.cdna.all.fa.gz', 'ref/Mus_musculus.GRCm38.ncrna.fa.gz'):
    for name, s in read_fasta(fn):
        tid = name.split()[0].split('.')[0]
        m = re.search(r'gene:(\S+)', name); gene_of[tid] = m.group(1).split('.')[0] if m else tid
        m = re.search(r'transcript_biotype:(\S+)', name); btype[tid] = m.group(1) if m else 'NA'
        seqs[tid] = s
TID = sorted(seqs)
tix = {t: i for i, t in enumerate(TID)}
print('tx', len(TID), round(time.time() - t0), flush=True)
var = {l: np.load(f'w/ref/var_{l}.npy') for l, _ in FOUNDERS}
# chrom lengths: max coordinate seen in GTF/variants
sel_codes, sel_tx, sel_pos = [], [], []
for ch in CHROMS:
    cnum = 20 if ch == 'X' else int(ch)
    tl = [t for t in TID if t in tx and tx[t]['chrom'] == ch and len(seqs[t]) == sum(e - s + 1 for s, e in tx[t]['exons'])]
    mx = max(max(e for s, e in tx[t]['exons']) for t in tl) + 10
    arr = np.zeros(mx + 2, np.uint8)
    for l, v in var.items():
        vv = v[(v[:, 0] == cnum) & (v[:, 1] < mx)]
        d = np.bincount(vv[:, 1], minlength=len(arr) + 1)[:len(arr) + 1].astype(np.int32)
        d -= np.bincount(np.minimum(vv[:, 2] + 1, len(arr)), minlength=len(arr) + 1)[:len(arr) + 1].astype(np.int32)
        arr |= (np.cumsum(d[:len(arr)]) > 0).astype(np.uint8) * np.uint8(FBIT[l])
    for t in tl:
        f, r, v = seqcodes(enc(seqs[t]))
        if not len(f):
            continue
        gc = tx_genome_coords(tx[t])
        m = arr[gc]
        win = np.zeros(len(f), np.uint8)
        for k in range(K):
            win |= m[k:k + len(f)]
        ok = np.nonzero(v & (win == 0))[0]
        if not len(ok):
            continue
        w = gc[ok] // SP
        _, first = np.unique(w, return_index=True)
        p = ok[np.sort(first)]
        sel_codes.append(f[p]); sel_tx.append(np.full(len(p), tix[t], np.int32)); sel_pos.append(p.astype(np.int32))
    print(ch, len(tl), round(time.time() - t0), flush=True)
B = np.unique(np.concatenate(sel_codes))
P26 = np.load('/mnt/project-files/revision/respacing/probes/probes_s26.npz')
A = np.array([seqcodes(enc(k))[0][0] for k in P26['kmer'].tolist()], np.uint64)
U = np.unique(np.concatenate([A, B]))
print('B', len(B), 'A', len(A), 'union', len(U), flush=True)
# compatibility: every sense occurrence of a union k-mer in any transcript
cp, ct, cpos = [], [], []
for t in TID:
    f, r, v = seqcodes(enc(seqs[t]))
    if not len(f):
        continue
    j = np.searchsorted(U, f); j[j == len(U)] = 0
    h = np.nonzero((U[j] == f) & v)[0]
    cp.append(j[h].astype(np.int32)); ct.append(np.full(len(h), tix[t], np.int32)); cpos.append(h.astype(np.int32))
cp, ct, cpos = (np.concatenate(x) for x in (cp, ct, cpos))
print('compat', len(cp), round(time.time() - t0), flush=True)
iA = np.searchsorted(U, A)
np.savez('design.npz', U=U, inB=np.isin(U, B), A_idx=iA, A_gene=P26['gene'], A_pos=P26['pos'], A_vf=P26['in_vf'], A_ref=P26['in_ref'],
         A_mask=P26['mask'], cp=cp, ct=ct, cpos=cpos, tid=np.array(TID), tgene=np.array([gene_of[t] for t in TID]),
         tlen=np.array([len(seqs[t]) for t in TID]), tbio=np.array([btype[t] for t in TID]),
         tkeep=np.array([bool(tx.get(t, {}).get('keep', False)) for t in TID]))
# antisense k-mers for KMC (read 1 of a reverse-stranded library)
def dec(c):
    out = []
    for x in c.tolist():
        s = []
        for i in range(K):
            s.append('ACGT'[(x >> (2 * (K - 1 - i))) & 3])
        out.append(''.join(s))
    return out
RC = str.maketrans('ACGT', 'TGCA')
with open('union_as.fa', 'w') as fo:
    for i, s in enumerate(dec(U)):
        fo.write(f'>{i}\n{s.translate(RC)[::-1]}\n')
print('done', round(time.time() - t0))
