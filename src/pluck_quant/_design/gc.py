import numpy as np
U = np.load('design.npz')['U']; x = U.copy(); gc = np.zeros(len(U), np.int64)
for i in range(25):
    b = x & np.uint64(3); gc += ((b == 1) | (b == 2)); x >>= np.uint64(2)
np.save('gc.npy', gc)
