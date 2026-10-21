#!/bin/bash
# Stage 4: Kmerator gene-level specific k-mers for the 16,885 core genes, then the tag list files for countTags.
# Needs km/genome_k31.jf (jellyfish --if index of the genome restricted to the cDNA+ncRNA k-mers; equivalence tested on chr19), km/DS, km/core_genes.txt.
set -uo pipefail
W=${W:-/home/user/fc}; cd $W/km; export PATH=$W/mm/envs/tools/bin:$PATH
s=$(date +%s)
[ -s KM_full/kmers.fa ] || [ -s share/tags.sha256 ] || kmerator -s core_genes.txt -g genome_k31.jf -d DS -S mouse -r 102 -k 31 -t ${KT:-3} -o KM_full -y > ../logs/km_full.log 2>&1
echo "kmerator_full $(( $(date +%s)-s )) s" >> ../logs/km_time.txt
SH=$W/km/share; mkdir -p $SH
{ [ -s $SH/tags.sha256 ] || python3 $W/scripts/mktags.py build KM_full/kmers.fa KM_full/contigs.fa $SH > ../logs/mktags_build.log 2>&1; } && \
python3 $W/scripts/mktags.py make $SH/contigs.fa.gz $SH $W/km/tags > ../logs/mktags_make.log 2>&1 && \
{ python3 -c "import random;random.seed(1);print(''.join('@r%d\n%s\n+\n%s\n'%(i,''.join(random.choice('ACGT') for _ in range(50)),'I'*50) for i in range(200)),end='')" > tags/dummy.fq; $W/src/countTags/bin/countTags -k 31 -i tags/tags.u.raw.gz --summary tags/dummy.sum tags/dummy.fq > tags/dummy.tsv 2>/dev/null && $W/mm/envs/tools/bin/python $W/scripts/fcsum.py ctperm tags/dummy.tsv tags/tags.u.raw.gz tags > ../logs/ctperm.log 2>&1 && rm tags/dummy.tsv && echo S4_KM_DONE; } || echo S4_KM_FAILED
