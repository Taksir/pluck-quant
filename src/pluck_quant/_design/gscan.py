# Genome occurrences (both strands) of each union probe, numba rolling 2-bit code + binary search on sorted canonical codes.
import gzip, numpy as np, numba, time
K = 25
D = np.load('design.npz'); U = D['U']
def rc(c):
    r = np.zeros_like(c)
    x = c.copy()
    for i in range(K):
        r = (r << np.uint64(2)) | (np.uint64(3) - (x & np.uint64(3))); x >>= np.uint64(2)
    return r
can = np.minimum(U, rc(U)); order = np.argsort(can); cs = can[order]
@numba.njit(cache=True)
def scan(b, cs, cnt):
    mask = (np.uint64(1) << np.uint64(2 * K)) - np.uint64(1)
    f = np.uint64(0); r = np.uint64(0); valid = 0
    n = len(cs)
    for i in range(len(b)):
        x = b[i]
        if x > 3:
            valid = 0; f = np.uint64(0); r = np.uint64(0); continue
        f = ((f << np.uint64(2)) | np.uint64(x)) & mask
        r = (r >> np.uint64(2)) | (np.uint64(3 - x) << np.uint64(2 * (K - 1)))
        valid += 1
        if valid >= K:
            c = f if f < r else r
            lo, hi = 0, n
            while lo < hi:
                mid = (lo + hi) >> 1
                if cs[mid] < c: lo = mid + 1
                else: hi = mid
            if lo < n and cs[lo] == c:
                cnt[lo] += 1
LUT = np.full(256, 4, np.uint8)
for i, ch in enumerate(b'ACGT'): LUT[ch] = i; LUT[ch + 32] = i
cnt = np.zeros(len(cs), np.int64)
t0 = time.time(); name = None; parts = []
def flush():
    if parts:
        b = LUT[np.frombuffer(b''.join(parts), np.uint8)]
        scan(b, cs, cnt); print(name, len(b), round(time.time() - t0), flush=True)
with gzip.open('ref/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz', 'rb') as fh:
    for line in fh:
        if line[:1] == b'>':
            flush(); name = line[1:].split()[0].decode(); parts = []
        else:
            parts.append(line.rstrip())
flush()
g = np.zeros(len(U), np.int64); g[order] = cnt
np.save('genome_count.npy', g)
print('GSCANDONE', (g == 0).mean(), (g == 1).mean(), (g > 1).mean())
