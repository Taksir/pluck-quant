#!/usr/bin/env python3
"""mktags.py build KMERS.fa CONTIGS.fa SHAREDIR   |   mktags.py make CONTIGS.fa SHAREDIR OUTDIR
(build: derive the exact pair set from Kmerator's kmers.fa, write the shareable contigs.fa.gz + patch_add/patch_del (pairs where kmers.fa and the contig windows differ, ~0.1%) + tags.sha256;
 make: regenerate the pair set from contigs.fa.gz + patches, verify tags.sha256, write the countTags files.)
Kmerator writes, per run, kmers.fa (every gene-specific 31-mer) and contigs.fa (maximal runs of consecutive specific 31-mers, >GENE:ENST.contig_N).
The tag list is regenerated from contigs.fa (small, shareable): all 31-mer windows of every contig, canonical (min of k-mer and reverse complement),
deduplicated per gene.  Outputs in OUTDIR:
  tags.u.raw.gz      unique canonical k-mers across all genes, sorted by 2-bit code, one per line  (countTags input -i)
  tag_uidx.npy       int32: for each (gene, k-mer) pair in gene order, its line index in tags.u.raw.gz
  tag_offsets.npy    int64: pair offset of each gene (len n_genes+1);  tag_genes.txt: gene ids in that order
  tags.sha256        sha256 of the text 'GENE<TAB>KMER<LF>' of all pairs sorted by (gene id, 2-bit code of the canonical k-mer)
  tags_stats.json
If KMERS.fa is given, the pair set built from it must equal the contigs-derived set (checked, exit 1 otherwise)."""
import sys, gzip, hashlib, json, numpy as np
K = 31
mode = sys.argv[1]
LUT = np.full(256, 255, np.uint8); LUT[np.frombuffer(b'ACGT', np.uint8)] = np.arange(4)
def read_contigs(path):
    genes, seqs = [], []
    with open(path) as f:
        for line in f:
            if line[0] == '>': genes.append(line[1:19])
            else: seqs.append(line.strip())
    assert len(genes) == len(seqs)
    return genes, seqs
def canon_pairs_from_contigs(genes, seqs):
    gid = {g: i for i, g in enumerate(sorted(set(genes)))}
    lens = np.array([len(s) for s in seqs], np.int64)
    b = LUT[np.frombuffer(''.join(seqs).encode(), np.uint8)]
    assert (b < 4).all()
    n = len(b); starts = np.concatenate([[0], np.cumsum(lens)[:-1]])
    cid = np.repeat(np.arange(len(seqs), dtype=np.int32), lens); pos = np.arange(n, dtype=np.int64) - starts[cid]
    ok = (pos <= lens[cid] - K) & (np.arange(n) <= n - K); del pos
    p = np.nonzero(ok)[0]; del ok
    can = np.empty(len(p), np.uint64); CH = 4_000_000
    for a in range(0, len(p), CH):
        pp = p[a:a + CH]; f = np.zeros(len(pp), np.uint64); r = np.zeros(len(pp), np.uint64)
        for j in range(K):
            x = b[pp + j].astype(np.uint64)
            f = (f << np.uint64(2)) | x
            r |= (np.uint64(3) - x) << np.uint64(2 * j)
        can[a:a + CH] = np.minimum(f, r)
    g = np.array([gid[x] for x in genes], np.int32)[cid[p]]; del cid, p
    o = np.lexsort((can, g)); g, can = g[o], can[o]
    keep = np.ones(len(g), bool); keep[1:] = (g[1:] != g[:-1]) | (can[1:] != can[:-1])
    return sorted(gid), g[keep], can[keep]
def decode(codes):
    c = codes.astype(np.uint64); m = np.empty((len(c), K), np.uint8)
    for j in range(K - 1, -1, -1): m[:, j] = np.frombuffer(b'ACGT', np.uint8)[(c & np.uint64(3)).astype(np.intp)]; c = c >> np.uint64(2)
    return m

def txt(g, can, ga):
    d = decode(can); m = np.empty((len(d), 18 + 1 + K + 1), np.uint8)
    m[:, :18] = np.frombuffer(ga[g].tobytes(), np.uint8).reshape(-1, 18); m[:, 18] = 9; m[:, 19:19 + K] = d; m[:, -1] = 10
    return m.tobytes()
def encode_kmers(ks):
    b = LUT[np.frombuffer(''.join(ks).encode(), np.uint8)].reshape(-1, K).astype(np.uint64); f = np.zeros(len(b), np.uint64); r = np.zeros(len(b), np.uint64)
    for j in range(K): f = (f << np.uint64(2)) | b[:, j]; r |= (np.uint64(3) - b[:, j]) << np.uint64(2 * j)
    return np.minimum(f, r)
def read_kmers(path, gl):
    gid = {x: i for i, x in enumerate(gl)}; gs, cs = [], []; kg, ks = [], []
    def flush():
        gs.append(np.array([gid[x] for x in kg], np.int32)); cs.append(encode_kmers(ks)); kg.clear(); ks.clear()
    with open(path) as f:
        for line in f:
            if line[0] == '>': kg.append(line[1:19])
            else:
                ks.append(line.strip())
                if len(ks) >= 2_000_000: flush()
    if ks: flush()
    g = np.concatenate(gs); c = np.concatenate(cs); del gs, cs
    o = np.lexsort((c, g)); g, c = g[o], c[o]
    keep = np.ones(len(g), bool); keep[1:] = (g[1:] != g[:-1]) | (c[1:] != c[:-1])
    return gl, g[keep], c[keep]
def read_patch(path, gl):
    gid = {x: i for i, x in enumerate(gl)}; gs, ks = [], []
    for line in gzip.open(path, 'rt'):
        a, b = line.rstrip('\n').split('\t'); gs.append(gid[a]); ks.append(b)
    return np.array(gs, np.int32), (encode_kmers(ks) if ks else np.zeros(0, np.uint64))
def minus(g, c, g2, c2):   # pairs of (g,c) not in (g2,c2); both inputs duplicate-free
    n = len(g); G = np.concatenate([g, g2]); C = np.concatenate([c, c2])
    o = np.lexsort((C, G)); Gs, Cs = G[o], C[o]; dup = np.zeros(len(o), bool)
    same = (Gs[1:] == Gs[:-1]) & (Cs[1:] == Cs[:-1]); dup[1:] |= same; dup[:-1] |= same
    out = np.ones(n, bool); idx = o[dup]; out[idx[idx < n]] = False
    return g[out], c[out]
if mode == 'build':
    kmers, contigs, shared = sys.argv[2:5]
    genes, seqs = read_contigs(contigs); gl = sorted(set(genes)); gl2, g2, c2 = canon_pairs_from_contigs(genes, seqs); del genes, seqs
    assert gl == gl2
    gl, g, c = read_kmers(kmers, gl)
    o = np.lexsort((c2, g2)); g2, c2 = g2[o], c2[o]
    ga = np.array(gl, dtype='S18')
    add_g, add_c = minus(g, c, g2, c2); del_g, del_c = minus(g2, c2, g, c)
    import shutil
    with open(contigs, 'rb') as fi, gzip.open(f'{shared}/contigs.fa.gz', 'wb', compresslevel=9) as fo: shutil.copyfileobj(fi, fo)
    for name, gg, cc in (('patch_add', add_g, add_c), ('patch_del', del_g, del_c)):
        with gzip.open(f'{shared}/{name}.txt.gz', 'wb') as fo:
            if len(gg): fo.write(txt(gg, cc, ga))
    h = hashlib.sha256(); B = 2_000_000
    for i in range(0, len(g), B): h.update(txt(g[i:i + B], c[i:i + B], ga))
    open(f'{shared}/tags.sha256', 'w').write(h.hexdigest() + '  tags (sorted GENE<TAB>canonical k-mer lines)\n')
    json.dump({'n_pairs': int(len(g)), 'n_genes': len(gl), 'patch_add': int(len(add_g)), 'patch_del': int(len(del_g)), 'sha256': h.hexdigest()}, open(f'{shared}/tags_stats.json', 'w'))
    print('build', len(g), 'pairs; patch add/del', len(add_g), len(del_g), h.hexdigest())
else:
    contigs, shared, out = sys.argv[2:5]
    op = gzip.open if contigs.endswith('.gz') else open
    import io; tmp = io.StringIO(op(contigs, 'rt').read())
    genes, seqs = [], []
    for line in tmp:
        if line[0] == '>': genes.append(line[1:19])
        else: seqs.append(line.strip())
    gl2, g2, c2 = canon_pairs_from_contigs(genes, seqs)
    want = open(f'{shared}/tags.sha256').read().split()[0]
    pg, pc = read_patch(f'{shared}/patch_add.txt.gz', None) if False else (None, None)
    # gene list of the final set = genes of kmers.fa; patches may reference genes absent from the contigs
    allg = set(gl2)
    for name in ('patch_add', 'patch_del'):
        for line in gzip.open(f'{shared}/{name}.txt.gz', 'rt'): allg.add(line.split('\t')[0])
    gl = sorted(allg); gid = {x: i for i, x in enumerate(gl)}
    g = np.array([gid[x] for x in gl2], np.int32)[g2]; c = c2
    dg, dc = read_patch(f'{shared}/patch_del.txt.gz', gl); ag, ac = read_patch(f'{shared}/patch_add.txt.gz', gl)
    if len(dg): g, c = minus(g, c, dg, dc)
    g = np.concatenate([g, ag]); c = np.concatenate([c, ac]); o = np.lexsort((c, g)); g, c = g[o], c[o]
    keep = np.ones(len(g), bool); keep[1:] = (g[1:] != g[:-1]) | (c[1:] != c[:-1]); g, c = g[keep], c[keep]
    ga = np.array(gl, dtype='S18'); h = hashlib.sha256(); B = 2_000_000
    for i in range(0, len(g), B): h.update(txt(g[i:i + B], c[i:i + B], ga))
    if h.hexdigest() != want: print('HASH MISMATCH', h.hexdigest(), want); sys.exit(1)
    ukm = np.unique(c); uidx = np.searchsorted(ukm, c).astype(np.int32); off = np.searchsorted(g, np.arange(len(gl) + 1)).astype(np.int64)
    with gzip.open(f'{out}/tags.u.raw.gz', 'wb', compresslevel=4) as fo:
        for i in range(0, len(ukm), B):
            d = decode(ukm[i:i + B]); m = np.empty((len(d), K + 1), np.uint8); m[:, :K] = d; m[:, K] = 10; fo.write(m.tobytes())
    np.save(f'{out}/tag_uidx.npy', uidx); np.save(f'{out}/tag_offsets.npy', off); open(f'{out}/tag_genes.txt', 'w').write('\n'.join(gl) + '\n')
    open(f'{out}/tags.sha256', 'w').write(h.hexdigest() + '  tags\n')
    print('make OK', len(g), 'pairs', len(gl), 'genes', len(ukm), 'unique k-mers; sha256 matches', want)
