"""Stage 3: one record per (gene, candidate 25-mer) for every gene on chr1-19,X, all transcripts.
Candidate = canonical 25-mer in <= GMAX genes of the Ensembl 102 transcriptome and <= GMAX times in GRCm38.
Fields: gene, can, sense (k-mer as in this gene's mRNA), meta (exon-union coordinate), nT1/nTA containing tx,
const/any flags, margin flags (>=100 nt from both ends of some containing tx), founder-variant mask (OR over occurrences)."""
import sys, pickle, time, numpy as np
sys.argv = ['x']; import e9rs as E
K, MARG = 25, 100
t0 = time.time()
A = pickle.load(open('w/annot.pkl', 'rb')); seqs = pickle.load(open('w/seqs.pkl', 'rb'))
tx, ok, genes_all = A['tx'], A['ok'], A['genes_all']
gid = {g: i for i, g in enumerate(genes_all)}
Z = np.load('w/kmer_ngene.npz'); UC, NG = Z['code'], Z['ngene']
U = np.load('w/U.npy'); UG = np.load('w/U_genome.npy')
var = {l: np.load(f'ref/var_{l}.npy') for l, _ in E.FOUNDERS}
genes = {}
for t in ok:
    genes.setdefault(tx[t]['gene_id'], []).append(t)
bychrom = {}
for g, ts in genes.items():
    bychrom.setdefault(tx[ts[0]]['chrom'], []).append(g)
chrom_len = {}
for line in open('ref/dna.fa.fai') if False else []:
    pass
out = {k: [] for k in ('gene', 'can', 'sense', 'meta', 'nT1', 'nTA', 'totT1', 'totTA', 'mT1', 'mTA', 'mask', 'ngene', 'gc')}
for ch in E.CHROMS:
    cnum = 20 if ch == 'X' else int(ch)
    gl = bychrom.get(ch, [])
    mx = max(max(e for t in genes[g] for _, e in tx[t]['exons']) for g in gl) + 10
    arr = np.zeros(mx + 2, np.uint8)
    for l, v in var.items():
        vv = v[v[:, 0] == cnum]
        vv = vv[vv[:, 1] < mx]
        d = np.zeros(mx + 3, np.int32)
        np.add.at(d, vv[:, 1], 1); np.add.at(d, np.minimum(vv[:, 2] + 1, mx + 2), -1)
        arr |= (np.cumsum(d[:mx + 2]) > 0).astype(np.uint8) * np.uint8(E.FBIT[l])
    for g in gl:
        ts = genes[g]
        strand = tx[ts[0]]['strand']
        gcs = [E.tx_genome_coords(tx[t]) for t in ts]
        allb = np.unique(np.concatenate(gcs))
        if strand == '-':
            allb = allb[::-1]
        meta_of = lambda gpos: (np.searchsorted(allb, gpos) if strand == '+' else len(allb) - 1 - np.searchsorted(allb[::-1], gpos))
        cols = {k: [] for k in ('can', 'f', 'ti', 'meta', 'mask', 'marg')}
        for i, (t, gc) in enumerate(zip(ts, gcs)):
            s = seqs[t]; L = len(s)
            f, r, v = E.seqcodes(E.enc(s))
            if not len(f):
                continue
            can = np.minimum(f, r)
            j = np.minimum(np.searchsorted(U, can), len(U) - 1)
            keep = v & (U[j] == can)
            keep[keep] &= UG[j[keep]] <= 8
            p = np.nonzero(keep)[0]
            if not len(p):
                continue
            m = arr[gc]
            win = np.zeros(len(m) - K + 1, np.uint8)
            for k in range(K):
                win |= m[k:k + len(win)]
            cols['can'].append(can[p]); cols['f'].append(f[p]); cols['ti'].append(np.full(len(p), i, np.int32))
            cols['meta'].append(meta_of(gc[p]).astype(np.int32)); cols['mask'].append(win[p])
            cols['marg'].append((p >= MARG) & (p <= L - MARG - K))
        if not cols['can']:
            continue
        c = {k: np.concatenate(v) for k, v in cols.items()}
        isT1 = np.array([tx[t]['keep'] for t in ts])
        o = np.lexsort((c['ti'], c['can'])); c = {k: v[o] for k, v in c.items()}
        st = np.r_[0, np.nonzero(c['can'][1:] != c['can'][:-1])[0] + 1]
        # distinct transcripts per k-mer
        newt = np.r_[True, (c['can'][1:] != c['can'][:-1]) | (c['ti'][1:] != c['ti'][:-1])]
        t1 = isT1[c['ti']]
        nTA = np.add.reduceat(newt.astype(np.int32), st)
        nT1 = np.add.reduceat((newt & t1).astype(np.int32), st)
        mT1 = np.maximum.reduceat((c['marg'] & t1).astype(np.uint8), st)
        mTA = np.maximum.reduceat(c['marg'].astype(np.uint8), st)
        mask = np.bitwise_or.reduceat(c['mask'], st)
        meta = np.minimum.reduceat(c['meta'], st)
        can = c['can'][st]
        j = np.searchsorted(UC, can); jj = np.searchsorted(U, can)
        out['gene'].append(np.full(len(st), gid[g], np.int32)); out['can'].append(can); out['sense'].append(c['f'][st])
        out['meta'].append(meta); out['nT1'].append(nT1.astype(np.int16)); out['nTA'].append(nTA.astype(np.int16))
        out['totT1'].append(np.full(len(st), isT1.sum(), np.int16)); out['totTA'].append(np.full(len(st), len(ts), np.int16))
        out['mT1'].append(mT1.astype(bool)); out['mTA'].append(mTA.astype(bool)); out['mask'].append(mask)
        out['ngene'].append(NG[j]); out['gc'].append(UG[jj])
    del arr
    print(ch, len(gl), sum(len(x) for x in out['gene']), round(time.time() - t0), flush=True)
np.savez('w/records.npz', **{k: np.concatenate(v) for k, v in out.items()})
# per-gene info for the model
info = []
for g, ts in genes.items():
    t = tx[ts[0]]
    gc = np.concatenate([E.tx_genome_coords(tx[x]) for x in ts])
    info.append(dict(gid=gid[g], gene_id=g, gene_name=t['gene_name'], chrom=t['chrom'], strand=t['strand'],
                     mid=int((gc.min() + gc.max()) // 2), tx=ts, tx_len=[len(seqs[x]) for x in ts], tx_t1=[bool(tx[x]['keep']) for x in ts]))
pickle.dump(info, open('w/geneinfo.pkl', 'wb'))
print('S3DONE', round(time.time() - t0), flush=True)
