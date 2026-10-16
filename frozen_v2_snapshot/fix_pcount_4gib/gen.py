"""Synthetic 101 bp SE FASTQ chunk for counter tests: random bases with 2 real-panel 25-mers planted per read, ~1% reads with an N.
usage: gen.py panel.txt NREADS SEED > chunk.fq"""
import sys, numpy as np
pan = np.frombuffer(open(sys.argv[1], 'rb').read(), np.uint8).reshape(-1, 26)[:, :25]
n, rng = int(sys.argv[2]), np.random.default_rng(int(sys.argv[3]))
ACGT = np.frombuffer(b'ACGT', np.uint8); out = sys.stdout.buffer; B = 200000
for s in range(0, n, B):
    m = min(B, n - s)
    seq = ACGT[rng.integers(0, 4, (m, 101))]
    for off in (rng.integers(0, 26, m), rng.integers(51, 77, m)):
        k = pan[rng.integers(0, len(pan), m)]
        idx = off[:, None] + np.arange(25); seq[np.arange(m)[:, None], idx] = k
    nn = rng.random(m) < 0.01; seq[nn, rng.integers(0, 101, nn.sum())] = ord('N')
    q = rng.choice(np.frombuffer(b'FFFF:,', np.uint8), (m, 101))
    lines = []
    for i in range(m):
        lines.append(b'@SIM.%d %d length=101\n%s\n+\n%s\n' % (s + i + 1, s + i + 1, seq[i].tobytes(), q[i].tobytes()))
    out.write(b''.join(lines))
