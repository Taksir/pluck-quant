#!/bin/bash
cd "${FR:-$PWD}"
c(){ [ -s "$2" ] || (curl -fsSL --retry 5 --retry-delay 3 --retry-all-errors -o "$2.part" "$1" && mv "$2.part" "$2"); }
c https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz e4/ref/genes.gtf.gz
ln -sf /home/user/fc/frames/e4/ref/genes.gtf.gz e5/w/ref/genes.gtf.gz
c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45684/suppl/GSE45684_FounderStrainsamples_alignedtoStraingenomes_processedmatrix.txt.gz" e5/dl/GSE45684_FounderStrainsamples_alignedtoStraingenomes_processedmatrix.txt.gz
while read f; do
  g=${f%%_*}
  c "https://ftp.ncbi.nlm.nih.gov/geo/samples/${g:0:7}nnn/$g/suppl/$f" e4/truth/$f
done < truth_files.txt
echo DONE > dl.done
