#!/usr/bin/env python3
"""gfix.py LAYOUT DESIGN_DIR : make genome_count.npy independent of numpy's argsort tie order.
Background: the frozen gscan.py sorts the canonical 25-mer codes with np.argsort (unstable); two probes with the same canonical code (a probe and its reverse complement)
form a tie group, the whole genome count lands on ONE member and the other gets 0, and which member depends on the sort implementation (AVX512 nodes differ from the machine
that built the locked designs).  The locked design has the carrier recorded in the original genome_count.npy; code/carriers{50,101}.b64 hold, for every tie group (ascending canonical code),
the rank (by index) of that carrier among the group's members.  This script keeps every group's TOTAL (invariant) and puts it on that member.  Idempotent.
Exit 1 if the number of tie groups differs from the table."""
import sys, zlib, base64, os, numpy as np
lay, d = sys.argv[1:3]; here = os.path.dirname(os.path.abspath(__file__)); K = 25
U = np.load(f'{d}/design.npz')['U']; g = np.load(f'{d}/genome_count.npy')
rank = np.frombuffer(zlib.decompress(base64.b64decode(open(f'{here}/carriers{lay}.b64').read())), np.uint8)
def rc(c):
    r = np.zeros_like(c); x = c.copy()
    for i in range(K):
        r = (r << np.uint64(2)) | (np.uint64(3) - (x & np.uint64(3))); x >>= np.uint64(2)
    return r
can = np.minimum(U, rc(U)); u, inv, cnt = np.unique(can, return_inverse=True, return_counts=True)
tg = np.nonzero(cnt > 1)[0]                       # tie groups, ascending canonical code
if len(tg) != len(rank): print('GFIX FAIL: tie groups', len(tg), 'table', len(rank)); sys.exit(1)
tot = np.bincount(inv, weights=g.astype(np.float64)); assert tot.max() < 2**52
gid = np.full(len(u), -1, np.int64); gid[tg] = np.arange(len(tg))
mem = np.nonzero(gid[inv] >= 0)[0]                # members of tie groups, ascending index
o = np.argsort(inv[mem], kind='stable'); mem = mem[o]   # grouped, ascending index inside a group
grp = inv[mem]; start = np.r_[0, np.nonzero(grp[1:] != grp[:-1])[0] + 1]
pos = np.arange(len(mem)) - np.repeat(start, np.diff(np.r_[start, len(mem)]))
carrier = pos == rank[gid[grp]]
new = g.copy(); new[mem] = np.where(carrier, tot[grp], 0).astype(np.int64)
print(f'GFIX layout {lay}: {len(tg)} tie groups, {int((new != g).sum())} members changed')
if (new != g).any(): np.save(f'{d}/genome_count.npy', new)
