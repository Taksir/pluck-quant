"""Stage 1: annotation, transcript sequences, transcriptome-wide gene count per canonical 25-mer.
Outputs w/annot.pkl, w/kmer_ngene.npz (unique codes + ngene), w/shared_pairs.npz (code, gene for 2..GMAX genes), w/U.npy."""
import sys, re, pickle, time, numpy as np
sys.argv = ['x']; import e9rs as E
GMAX = 8
t0 = time.time()
tx = E.parse_gtf('ref/genes.gtf.gz')            # CHROMS only, exons, keep flag = TSL1 pc
# biotype/tsl not kept by parse_gtf beyond keep; fine
seqs, gene_of = {}, {}
for fn in ('ref/cdna.fa.gz', 'ref/ncrna.fa.gz'):
    for name, s in E.read_fasta(fn):
        tid = name.split()[0].split('.')[0]
        m = re.search(r'gene:(\S+)', name)
        gene_of[tid] = m.group(1).split('.')[0] if m else tid
        seqs[tid] = s
# universe transcripts: on CHROMS, sequence length matches exons
ok = {t for t, r in tx.items() if t in seqs and len(seqs[t]) == sum(e - s + 1 for s, e in r['exons'])}
genes_all = sorted(set(gene_of.values()))
gid = {g: i for i, g in enumerate(genes_all)}
print('transcripts', len(seqs), 'genes', len(genes_all), 'universe tx', len(ok), flush=True)
codes, gidx = [], []
per_gene = {}
for t, g in gene_of.items():
    per_gene.setdefault(g, []).append(t)
for g, ts in per_gene.items():
    cs = []
    for t in ts:
        f, r, v = E.seqcodes(E.enc(seqs[t]))
        cs.append(np.minimum(f, r)[v])
    u = np.unique(np.concatenate(cs))
    codes.append(u); gidx.append(np.full(len(u), gid[g], np.int32))
codes = np.concatenate(codes); gidx = np.concatenate(gidx)
o = np.argsort(codes, kind='stable'); codes = codes[o]; gidx = gidx[o]; del o
print('pairs', len(codes), round(time.time() - t0), flush=True)
st = np.r_[0, np.nonzero(codes[1:] != codes[:-1])[0] + 1]
uc = codes[st]; ng = np.diff(np.r_[st, len(codes)])
np.savez('w/kmer_ngene.npz', code=uc, ngene=np.minimum(ng, 255).astype(np.uint8))
shared = np.repeat((ng >= 2) & (ng <= GMAX), ng)
np.savez('w/shared_pairs.npz', code=codes[shared], gene=gidx[shared])
np.save('w/U.npy', uc[ng <= GMAX])
print('distinct', len(uc), 'ngene<=8', int((ng <= GMAX).sum()), 'ng1', int((ng == 1).sum()), 'ng2', int((ng == 2).sum()), flush=True)
pickle.dump(dict(tx=tx, ok=ok, gene_of=gene_of, genes_all=genes_all), open('w/annot.pkl', 'wb'))
pickle.dump(seqs, open('w/seqs.pkl', 'wb'))
print('S1DONE', round(time.time() - t0), flush=True)
