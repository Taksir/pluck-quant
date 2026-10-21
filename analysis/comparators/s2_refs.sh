#!/bin/bash
# Stage 2: Ensembl 102 + MGP v5 references and the founder transcriptomes (ref.fa = B6 only, pan.fa = B6 + founder copies, tx2gene.tsv).
set -uo pipefail
W=${W:-/home/user/fc}; cd $W; T=$W/mm/envs/tools/bin
ENS=https://ftp.ensembl.org/pub/release-102
dl(){ for i in 1 2 3 4 5; do [ -s $2 ] && return 0; $T/aria2c -q -x16 -s16 -k4M --max-tries=5 --allow-overwrite=true -d $(dirname $2) -o $(basename $2).part "$1" && { case $2 in *.gz) gzip -t $2.part || { rm -f $2.part; continue; };; esac; mv $2.part $2 && return 0; }; rm -f $2.part; sleep $((5*i)); done; echo "FAIL $2"; return 1; }
mkdir -p ref/vcf fa
dl $ENS/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz ref/cdna.fa.gz
dl $ENS/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz ref/ncrna.fa.gz
dl $ENS/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz ref/genes.gtf.gz
dl https://raw.githubusercontent.com/rqtl/qtl2data/main/CC/cc.zip ref/cc.zip
if [ ! -s fa/stats.json ]; then
  dl $ENS/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ref/dna.fa.gz
  M=https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs
  for s in A_J 129S1_SvImJ NOD_ShiLtJ NZO_HlLtJ CAST_EiJ PWK_PhJ WSB_EiJ; do for k in snps.dbSNP142.vcf.gz indels.dbSNP142.normed.vcf.gz; do
    dl $M/$s.mgp.v5.$k ref/vcf/$s.mgp.v5.$k & done; done; wait
  $T/python code/build_founders.py ref fa > logs/build_founders.log 2>&1 && echo founders_ok
fi
ls -la fa; cat fa/stats.json | head -c 400; echo
