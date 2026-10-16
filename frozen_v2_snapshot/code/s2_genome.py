"""Stage 2: occurrences (canonical, capped at 255) of each candidate 25-mer U in the GRCm38 primary assembly."""
import sys, re, time, numpy as np, multiprocessing as mp
sys.argv = ['x']; import e9rs as E
NW = 4; STEP = 10_000_000
def load_genome():
    data = open('ref/dna.fa', 'rb').read()
    st = [m.start() for m in re.finditer(rb'>', data)] + [len(data)]
    recs = []
    for i in range(len(st) - 1):
        rec = data[st[i]:st[i + 1]]; nl = rec.index(b'\n')
        recs.append(E.LUT[np.frombuffer(rec[nl + 1:].replace(b'\n', b''), np.uint8)])
    return recs
RECS = U = None
def work(wi):
    recs = RECS
    cnt = np.zeros(len(U), np.uint32); k = 0
    for b in recs:
        for s in range(0, len(b), STEP):
            k += 1
            if k % NW != wi: continue
            f, r, v = E.seqcodes(b[s:s + STEP + E.K - 1])
            c = np.minimum(f, r)[v]
            j = np.searchsorted(U, c); j[j == len(U)] = 0
            h = j[U[j] == c]
            if len(h) < 200000:
                np.add.at(cnt, h, 1)
            else:
                cnt += np.bincount(h, minlength=len(U)).astype(np.uint32)
    return cnt
if __name__ == '__main__':
    t0 = time.time()
    import os, subprocess
    if not os.path.exists('ref/dna.fa'):
        subprocess.run("gzip -dc ref/dna.fa.gz > ref/dna.fa", shell=True, check=True)
    U = np.load('w/U.npy'); RECS = load_genome()
    print('loaded', round(time.time() - t0), flush=True)
    with mp.get_context('fork').Pool(NW) as p:
        parts = p.map(work, range(NW))
    tot = np.zeros(len(parts[0]), np.int64)
    for x in parts: tot += x
    np.save('w/U_genome.npy', np.minimum(tot, 255).astype(np.uint8))
    print('S2DONE', round(time.time() - t0), np.bincount(np.minimum(tot, 10)).tolist(), flush=True)
