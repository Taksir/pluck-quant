#!/usr/bin/env python3
"""RNA-Skim reference prep (no variants). Usage: rs_prep.py REFDIR OUTDIR [CHROM|-] [protein_coding]
REFDIR has cdna.fa.gz and ncrna.fa.gz (Ensembl 102, GRCm38; symlinks fine).
Writes OUTDIR/ref.fa (one record per transcript, id = ENSMUST without version),
OUTDIR/tx2gene.tsv (tid<TAB>gene id, version stripped), OUTDIR/gene.fa (RNA-Skim 'specialized' FASTA:
'>G|T1|T2' then 'seq1|seq2'). Same transcript set and gene mapping as e9_salmon/build_founders.py's ref.fa / tx2gene.tsv
(cdna then ncrna, later record overwrites an earlier one with the same id)."""
import gzip, os, re, sys, hashlib

def read_fasta(fn):
    name, buf = None, []
    with gzip.open(fn, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if name is not None: yield name, ''.join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None: yield name, ''.join(buf)

ref, out = sys.argv[1], sys.argv[2]
only_chrom = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != '-' else None   # test mode only: keep transcripts whose header location is on this chromosome
pc_only = len(sys.argv) > 4 and sys.argv[4] == 'protein_coding'   # keep transcripts of genes with gene_biotype:protein_coding (RNA-Skim README's recommended set; ncRNA file has none)
os.makedirs(out, exist_ok=True)
seqs, gene_of = {}, {}
for fn in ('cdna.fa.gz', 'ncrna.fa.gz'):
    for name, s in read_fasta(os.path.join(ref, fn)):
        tid = name.split()[0].split('.')[0]
        m = re.search(r'gene:(\S+)', name)
        gene_of[tid] = m.group(1).split('.')[0] if m else tid
        if only_chrom and not re.search(r'(?:chromosome|primary_assembly|scaffold|toplevel):GRCm38:%s:' % re.escape(only_chrom), name): continue
        if pc_only and 'gene_biotype:protein_coding' not in name: continue
        seqs[tid] = s.upper()
by_gene = {}
for tid in seqs:
    by_gene.setdefault(gene_of[tid], []).append(tid)
with open(os.path.join(out, 'ref.fa'), 'w') as f, open(os.path.join(out, 'tx2gene.tsv'), 'w') as m:
    for tid, s in seqs.items():
        f.write(f'>{tid}\n{s}\n'); m.write(f'{tid}\t{gene_of[tid]}\n')
with open(os.path.join(out, 'gene.fa'), 'w') as g:
    for gid, tids in by_gene.items():
        g.write('>' + '|'.join([gid] + tids) + '\n' + '|'.join(seqs[t] for t in tids) + '\n')
short = sum(1 for s in seqs.values() if len(s) < 60)
print(f'transcripts={len(seqs)} genes={len(by_gene)} shorter_than_60={short}')
