#!/usr/bin/env python3
"""Gene-level sums for the final comparison. All outputs: TSV.gz with columns gene, value (plus extras where stated).
  fcsum.py kallisto  RUN_DIR  tx2gene.tsv  OUT.tsv.gz      sum of est_counts per gene (abundance.tsv); run_info.json is copied by the runner
  fcsum.py salmon    quant.sf tx2gene.tsv  OUT.tsv.gz      sum of NumReads per gene
  fcsum.py fcounts   fc.txt   OUT.tsv.gz                   featureCounts gene counts (column 7)
  fcsum.py counttags CT.tsv TAGDIR OUT.tsv.gz  |  ctperm CT.tsv tags.u.raw.gz TAGDIR  per gene: mean count over its tags, sum, n tags, n tags with count>0
tx2gene.tsv: transcript_id<TAB>gene_id (founder copies like ENSMUST..._A map to the same gene)."""
import sys, gzip, collections
def t2g(path): return dict(l.rstrip('\n').split('\t')[:2] for l in open(path))
def write(path, header, rows):
    with gzip.open(path, 'wt') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
cmd = sys.argv[1]
if cmd in ('kallisto', 'salmon'):
    src, tx, out = sys.argv[2:5]; m = t2g(tx); s = collections.defaultdict(float); tp = collections.defaultdict(float); miss = 0
    path = f'{src}/abundance.tsv' if cmd == 'kallisto' else src
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as f:
        h = f.readline().rstrip('\n').split('\t'); ci = h.index('est_counts' if cmd == 'kallisto' else 'NumReads'); ti = h.index('tpm' if cmd == 'kallisto' else 'TPM')
        for line in f:
            x = line.rstrip('\n').split('\t'); g = m.get(x[0])
            if g is None: miss += 1; continue
            s[g] += float(x[ci]); tp[g] += float(x[ti])
    write(out, ['gene', 'est_counts', 'tpm'], [(g, '%.4f' % v, '%.4f' % tp[g]) for g, v in sorted(s.items())]); print(f'{cmd}: {len(s)} genes, {miss} unmapped transcript ids')
elif cmd == 'fcounts':
    src, out = sys.argv[2:4]; rows = []
    for line in open(src):
        if line[0] == '#' or line.startswith('Geneid'): continue
        x = line.rstrip('\n').split('\t'); rows.append((x[0], x[-1]))
    write(out, ['gene', 'count'], rows); print(f'fcounts: {len(rows)} genes')
elif cmd == 'ctperm':
    # fcsum.py ctperm CT.tsv(any countTags output on tags.u.raw.gz) tags.u.raw.gz OUTDIR : row order of countTags output -> line index in tags.u.raw.gz (countTags prints in a fixed hash order, and in its own canonical orientation)
    import numpy as np, subprocess, hashlib, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    ct, raw, outd = sys.argv[2:5]
    LUT = np.full(256, 255, np.uint8); LUT[np.frombuffer(b'ACGT', np.uint8)] = np.arange(4); K = 31
    def enc(b):
        b = LUT[b].reshape(-1, K).astype(np.uint64); f = np.zeros(len(b), np.uint64); r = np.zeros(len(b), np.uint64)
        for j in range(K): f = (f << np.uint64(2)) | b[:, j]; r |= (np.uint64(3) - b[:, j]) << np.uint64(2 * j)
        return np.minimum(f, r)
    def codes(cmd):
        out = []; p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, executable='/bin/bash')
        while True:
            buf = p.stdout.read(32 * 4000000)
            if not buf: break
            rest = b''
            if len(buf) % 32:
                more = p.stdout.read(32 - len(buf) % 32); buf += more
            a = np.frombuffer(buf, np.uint8).reshape(-1, 32)
            assert (a[:, 31] == 10).all(); out.append(enc(a[:, :31].copy().reshape(-1)))
        p.wait(); assert p.returncode == 0; return np.concatenate(out)
    u = codes(f'zcat {raw}')                      # sorted canonical codes, one per line
    assert (u[1:] > u[:-1]).all()
    c = codes(f'tail -n +2 {ct} | cut -f1')
    perm = np.searchsorted(u, c).astype(np.int32)
    assert len(c) == len(u) and (u[perm] == c).all() and len(np.unique(perm)) == len(u)
    os.makedirs(outd, exist_ok=True); np.save(f'{outd}/ct_perm.npy', perm)
    h = subprocess.run(f'tail -n +2 {ct} | cut -f1 | sha256sum', shell=True, capture_output=True, text=True, executable='/bin/bash').stdout.split()[0]
    open(f'{outd}/ct_order.sha256', 'w').write(h + '\n'); print('ctperm ok', len(u), h)
elif cmd == 'counttags':
    # fcsum.py counttags CT.tsv TAGDIR OUT.tsv.gz : per gene mean count over its specific tags (positional; order checked against ct_order.sha256)
    import numpy as np, subprocess, os
    src, td, out = sys.argv[2:5]
    perm = np.load(f'{td}/ct_perm.npy'); uidx = np.load(f'{td}/tag_uidx.npy'); off = np.load(f'{td}/tag_offsets.npy'); genes = open(f'{td}/tag_genes.txt').read().split()
    h = subprocess.run(f'tail -n +2 {src} | cut -f1 | sha256sum', shell=True, capture_output=True, text=True, executable='/bin/bash').stdout.split()[0]
    if not os.environ.get('FC_SKIP_ORDER_CHECK'): assert h == open(f'{td}/ct_order.sha256').read().strip(), 'countTags row order changed'
    cnt = np.fromstring(subprocess.run(f'tail -n +2 {src} | cut -f2', shell=True, capture_output=True, executable='/bin/bash').stdout, dtype=np.int64, sep='\n') if False else np.array(subprocess.run(f'tail -n +2 {src} | cut -f2', shell=True, capture_output=True, executable='/bin/bash').stdout.split(), dtype=np.int64)
    assert len(cnt) == len(perm)
    cu = np.zeros(len(perm), np.int64); cu[perm] = cnt
    v = cu[uidx]; n = np.diff(off); cs = np.concatenate([[0], np.cumsum(v)]); sm = cs[off[1:]] - cs[off[:-1]]
    pc = np.concatenate([[0], np.cumsum(v > 0)]); npos = pc[off[1:]] - pc[off[:-1]]
    rows = [(genes[i], '%.6f' % (sm[i] / n[i]), int(sm[i]), int(n[i]), int(npos[i])) for i in range(len(genes)) if n[i] > 0]
    write(out, ['gene', 'mean', 'sum', 'n_tags', 'n_pos'], rows); print(f'counttags: {len(rows)} genes, {int((cnt>0).sum())} tags with count>0')
