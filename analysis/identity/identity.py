"""Sample identity check: diagnostic k-mer counts from the panel pass (both orientations summed), then the score step of
dk_gw.py (genome-wide copy of dk.py). Writes P/identity.tsv: keep = labelled strain is the top fit and fit >= 0.90."""
import sys, os, json, subprocess, numpy as np
# Folder constants below point to the original working folder; set them for your copy.
W = '/home/user/e4'; P = '/mnt/project-files/revision/e4_kidney'
mp = np.load(f'{W}/map.npz'); n = int(mp['nsites']); f, r = mp['idk_f'], mp['idk_r']
srrs = [l.split('\t')[0] for l in open(f'{P}/prereg/samples.tsv').read().splitlines()[1:]]
srrs = [s for s in srrs if os.path.exists(f'{P}/done/{s}')]
os.makedirs(f'{W}/dk/counts', exist_ok=True)
for s in srrs:
    c = np.load(f'{P}/counts/{s}.npz')['c'].astype(np.int64)
    reads = json.load(open(f'{P}/salmon/{s}_ref.meta_info.json'))['num_processed']
    tot = c[f] + c[r]
    np.savez(f'{W}/dk/counts/{s}.npz', ref=tot[:n], alt=tot[n:], reads=reads, read_len=50)
subprocess.run([sys.executable, f'{W}/code/dk_gw.py', 'score', f'{W}/dk'] + srrs, check=True, stdout=open(f'{W}/logs/dkscore.log', 'w'))
S = json.load(open(f'{W}/dk/summary.json'))
with open(f'{P}/identity.tsv', 'w') as o:
    o.write('srr\tstrain\tgenes_called\tconsistent\tlabelled_fit\trank_of_labelled\ttop_strain\ttop_fit\tsecond_strain\tsecond_fit\tkeep\n')
    for s in srrs:
        x = S[s]; top = x['identity_top5']; lab = dict(top).get(x['strain'])
        rank = x['identity_rank_of_true']
        if lab is None:  # labelled strain not in top 5
            lab = float('nan')
        keep = rank == 1 and top[0][1] >= 0.90
        o.write(f"{s}\t{x['strain']}\t{x['genes_called']}\t{x['consistent']}\t{lab}\t{rank}\t{top[0][0]}\t{top[0][1]}\t{top[1][0]}\t{top[1][1]}\t{keep}\n")
os.system(f'cp {W}/dk/summary.json {P}/identity_summary.json')
print(open(f'{P}/identity.tsv').read())
