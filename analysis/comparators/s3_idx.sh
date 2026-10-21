#!/bin/bash
# Stage 3: indexes (kallisto ref/pan, Salmon ref/pan, STAR sparse-D). STAR runs alone (RAM). Timings -> logs/idx_timing.tsv. Resumable.
set -uo pipefail
W=${W:-/home/user/fc}; cd $W; T=$W/mm/envs/tools/bin; SAL=$W/salmon-1.10.0/bin/salmon
tm(){ local n=$1; shift; local a=$(date +%s.%N); "$@"; local rc=$?; printf "%s\t%.0f\t%s\n" $n $(echo "$(date +%s.%N) - $a" | bc) $rc >> logs/idx_timing.tsv; return $rc; }
mkdir -p idx
[ -s idx/kallisto_ref.idx ] || tm kallisto_ref $T/kallisto index -k 31 -t 4 -i idx/kallisto_ref.idx fa/ref.fa > logs/kidx_ref.log 2>&1
[ -s idx/kallisto_pan.idx ] || tm kallisto_pan $T/kallisto index -k 31 -t 4 -i idx/kallisto_pan.idx fa/pan.fa > logs/kidx_pan.log 2>&1
[ -s idx/salmon_ref/versionInfo.json ] || tm salmon_ref $SAL index -t fa/ref.fa -i idx/salmon_ref -k 31 -p 4 > logs/sidx_ref.log 2>&1
[ -s idx/salmon_pan/versionInfo.json ] || tm salmon_pan $SAL index -t fa/pan.fa -i idx/salmon_pan -k 31 -p 4 > logs/sidx_pan.log 2>&1
if [ ! -s idx/star_sparse/SA ]; then
  [ -s ref/dna.fa ] || gzip -dc ref/dna.fa.gz > ref/dna.fa
  [ -s ref/genes.gtf ] || gzip -dc ref/genes.gtf.gz > ref/genes.gtf
  mkdir -p idx/star_sparse
  tm star_sparse $T/STAR --runMode genomeGenerate --runThreadN 4 --genomeDir idx/star_sparse --genomeFastaFiles ref/dna.fa --sjdbGTFfile ref/genes.gtf \
     --sjdbOverhang 100 --genomeSAsparseD ${STAR_D:-10} --limitGenomeGenerateRAM 13000000000 --outFileNamePrefix logs/star_gen_ > logs/star_gen.out 2>&1
  rm -rf idx/star_sparse/_STARtmp
fi
du -sm idx/* | tee logs/idx_sizes.txt
[ -s idx/star_sparse/SA ] && [ -s idx/kallisto_pan.idx ] && [ -s idx/salmon_pan/versionInfo.json ] && echo S3_IDX_DONE || echo S3_IDX_FAILED
