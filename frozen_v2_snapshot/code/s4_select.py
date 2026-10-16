"""Stage 4: probe sets. Each set = eligibility rule on stage-3 records + 26 nt spacing on the gene's exon-union
coordinate. Baseline probes (respacing/probes/probes_s26.npz, the current 26 nt set) are placed first wherever they
are eligible; relaxed candidates fill the gaps >= 26 nt from every selected probe. ref = no variant filter, vf = no
founder SNP/indel at any annotated copy. Stage 5: probe x transcript compatibility (read starts covering the probe)."""
import sys, json, pickle, time, numpy as np
from numba import njit
sys.argv = ['x']; import e9rs as E
S, R, K = 26, 50, 25
t0 = time.time()
Rc = dict(np.load('w/records.npz'))
n = len(Rc['can']); print('records', n, flush=True)
info = pickle.load(open('w/geneinfo.pkl', 'rb'))
A = pickle.load(open('w/annot.pkl', 'rb')); genes_all = A['genes_all']
# ---- group stats over all universe genes containing the same canonical k-mer
o = np.argsort(Rc['can'], kind='stable')
cs = Rc['can'][o]
st = np.r_[0, np.nonzero(cs[1:] != cs[:-1])[0] + 1]; ln = np.diff(np.r_[st, n])
def bcast(vals):
    out = np.empty(n, vals.dtype); out[o] = np.repeat(vals, ln); return out
complete = bcast((ln == Rc['ngene'][o][st]))
maskor = bcast(np.bitwise_or.reduceat(Rc['mask'][o], st))
anyT1 = bcast(np.minimum.reduceat((Rc['nT1'][o] > 0).astype(np.uint8), st)).astype(bool)
constT1 = bcast(np.minimum.reduceat(((Rc['nT1'] == Rc['totT1']) & (Rc['totT1'] > 0)).astype(np.uint8)[o], st)).astype(bool)
constTA = bcast(np.minimum.reduceat((Rc['nTA'] == Rc['totTA']).astype(np.uint8)[o], st)).astype(bool)
del o, cs
# ---- baseline (current 26 nt set): match by canonical code among single-gene records
P = np.load('/mnt/project-files/revision/respacing/probes/probes_s26.npz')
G0 = json.load(open('/mnt/project-files/revision/fast_competitors/e9_prep_cache/ref/genes.json'))
def codes(kmers):
    b = E.enc(''.join(kmers)).reshape(-1, K)
    f = np.zeros(len(b), np.uint64)
    for i in range(K): f = (f << np.uint64(2)) | (b[:, i] & 3).astype(np.uint64)
    r = np.zeros(len(b), np.uint64)
    for i in range(K): r |= (np.uint64(3) - (b[:, i] & 3).astype(np.uint64)) << np.uint64(2 * i)
    return f, np.minimum(f, r)
pf, pc = codes(P['kmer'].tolist())
gidmap = {g: i for i, g in enumerate(genes_all)}
pg = np.array([gidmap[G0[g]['gene_id']] for g in P['gene']])
single = Rc['ngene'] == 1
base = {}
for v in ('ref', 'vf'):
    sel = P['in_' + v]
    key = np.unique(pc[sel])
    j = np.minimum(np.searchsorted(key, Rc['can']), len(key) - 1)
    base[v] = single & (key[j] == Rc['can'])
    print('baseline', v, int(sel.sum()), 'matched records', int(base[v].sum()), flush=True)
# ---- eligibility
isT1gene = Rc['totT1'] > 0
spec = {'S1': (Rc['ngene'] == 1) & (Rc['gc'] <= 1), 'S2': (Rc['ngene'] <= 2) & (Rc['gc'] <= 2), 'S8': complete}
T1 = complete & isT1gene & anyT1 & Rc['mT1']
TA = complete & Rc['mTA']
SETS = {
    'base':    None,
    'T1_S2_C': T1 & spec['S2'] & constT1,
    'T1_S1_P': T1 & spec['S1'],
    'T1_S2_P': T1 & spec['S2'],
    'TA_S1_C': TA & spec['S1'] & constTA,
    'TA_S1_P': TA & spec['S1'],
    'TA_S2_P': TA & spec['S2'],
    'TA_S8_P': TA & spec['S8'],
}
vfok = maskor == 0
# ---- spacing
order = np.lexsort((Rc['meta'], Rc['gene']))
gs = Rc['gene'][order]; meta = Rc['meta'][order].astype(np.int64)
gb = np.r_[0, np.nonzero(gs[1:] != gs[:-1])[0] + 1, len(gs)]
@njit(cache=True)
def space(gb, meta, B, El, S):
    sel = np.zeros(len(meta), np.bool_)
    s1 = np.empty(len(meta), np.int64)
    for gi in range(len(gb) - 1):
        a, b = gb[gi], gb[gi + 1]
        last = -10 ** 9; ns = 0
        for i in range(a, b):
            if B[i] and El[i] and meta[i] - last >= S:
                sel[i] = True; last = meta[i]; s1[ns] = meta[i]; ns += 1
        lastany = -10 ** 9; p = 0
        for i in range(a, b):
            if sel[i]:
                if meta[i] > lastany: lastany = meta[i]
                continue
            if not El[i]: continue
            while p < ns and s1[p] < meta[i]: p += 1
            nxt = s1[p] if p < ns else 10 ** 12
            if meta[i] - lastany >= S and nxt - meta[i] >= S:
                sel[i] = True; lastany = meta[i]
    return sel
chosen = {}
for name, El in SETS.items():
    for v in ('ref', 'vf'):
        if El is None:
            e = base[v]
        else:
            e = El & (vfok if v == 'vf' else True)
        sel = space(gb, meta, base[v][order], e[order], S)
        idx = order[sel]
        chosen[(name, v)] = idx
        del sel
        print(name, v, 'records', len(idx), 'genes', len(np.unique(Rc['gene'][idx])), round(time.time() - t0), flush=True)
# ---- probe union (sense strings), per-set membership
allidx = np.unique(np.concatenate(list(chosen.values())))
usense = np.unique(Rc['sense'][allidx])
print('union probes', len(usense), flush=True)
memb = {}
for (name, v), idx in chosen.items():
    memb[f'{name}__{v}'] = np.unique(np.searchsorted(usense, Rc['sense'][idx])).astype(np.int32)
# record-level detail for the union (for the EM thread): gene, ngene, gc, mask of every record whose sense is a probe
# ---- stage 5: compatibility over every universe transcript (same orientation only)
rs = np.searchsorted(usense, Rc['sense'][allidx])
pmask = np.zeros(len(usense), np.uint8); np.bitwise_or.at(pmask, rs, maskor[allidx])
png = np.zeros(len(usense), np.uint8); np.maximum.at(png, rs, Rc['ngene'][allidx])
pgc = np.zeros(len(usense), np.uint8); np.maximum.at(pgc, rs, Rc['gc'][allidx])
del Rc, maskor, complete, anyT1, constT1, constTA, order, gs, meta, base, SETS, spec, T1, TA, vfok, single
import gc; gc.collect()
seqs = pickle.load(open('w/seqs.pkl', 'rb'))
txl, txg, txlen, txt1 = [], [], [], []
for g in info:
    for t, L, t1 in zip(g['tx'], g['tx_len'], g['tx_t1']):
        txl.append(t); txg.append(g['gid']); txlen.append(L); txt1.append(t1)
oc_p, oc_t, oc_n = [], [], []
for ti, t in enumerate(txl):
    f, r, v = E.seqcodes(E.enc(seqs[t]))
    if not len(f): continue
    j = np.minimum(np.searchsorted(usense, f), len(usense) - 1)
    h = np.nonzero(v & (usense[j] == f))[0]
    if not len(h): continue
    L = len(seqs[t])
    nst = np.maximum(0, np.minimum(h, L - R) - np.maximum(0, h - (R - K)) + 1)
    oc_p.append(j[h]); oc_t.append(np.full(len(h), ti, np.int32)); oc_n.append(nst)
oc_p = np.concatenate(oc_p).astype(np.int32); oc_t = np.concatenate(oc_t); oc_n = np.concatenate(oc_n).astype(np.int32)
# merge repeated occurrences of a probe in one transcript
key = oc_p.astype(np.int64) * (len(txl) + 1) + oc_t
uk, inv = np.unique(key, return_inverse=True)
A_n = np.bincount(inv, weights=oc_n).astype(np.int32)
A_p = (uk // (len(txl) + 1)).astype(np.int32); A_t = (uk % (len(txl) + 1)).astype(np.int32)
# per-probe attributes (from the record(s) with that sense): founder mask OR, ngene, genome count
def tostr(c):
    b = np.zeros((len(c), K), np.uint8); x = c.copy()
    for i in range(K - 1, -1, -1):
        b[:, i] = np.frombuffer(b'ACGT', np.uint8)[(x & np.uint64(3)).astype(np.int64)]; x >>= np.uint64(2)
    return b
sb = tostr(usense)
anti = np.frombuffer(bytes(sb[:, ::-1].tobytes()).translate(bytes.maketrans(b'ACGT', b'TGCA')), np.uint8).reshape(-1, K)
with open('w/panel_antisense.txt', 'wb') as fh:
    fh.write(b'\n'.join(bytes(r) for r in anti) + b'\n')
np.savez_compressed('w/design.npz', sense_code=usense, sense=np.array([bytes(r).decode() for r in sb]),
                    mask=pmask, ngene=png, genome_count=pgc,
                    A_probe=A_p, A_tx=A_t, A_nstarts=A_n,
                    tx_id=np.array(txl), tx_gene=np.array(txg, np.int32), tx_len=np.array(txlen, np.int32), tx_tsl1pc=np.array(txt1),
                    gene_id=np.array(genes_all), **{'set_' + k: v for k, v in memb.items()})
print('S4DONE', len(usense), len(A_p), round(time.time() - t0), flush=True)
