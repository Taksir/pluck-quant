#!/bin/bash
# One accuracy-leg sample: download -> kallisto ref/pan -> countTags -> STAR+featureCounts -> seqtk 25%/10% -> Salmon ref/pan on each -> gene sums.
# usage: fc_sample.sh COHORT SRR        (rows come from $W/samples_fc.tsv)   env: W (work dir, default: the original working folder), THREADS (default 4)
# Local outputs: $W/out/COHORT/SRR/  (small files only; FASTQ/BAM/subsamples are deleted at the end). Idempotent: skips a sample with out/.../DONE.
set -uo pipefail
W=${W:-/home/user/fc}; TH=${THREADS:-4}; T=$W/mm/envs/tools/bin; SAL=$W/salmon-1.10.0/bin/salmon; CT=$W/src/countTags/bin/countTags
S=${S:-$W/scripts}; COH=$1; SRR=$2
row=$(awk -F'\t' -v c=$COH -v s=$SRR '$1==c && $2==s' $W/samples_fc.tsv); [ -n "$row" ] || { echo "no such sample $COH $SRR"; exit 2; }
IFS=$'\t' read -r _ _ title strain nreads ftp strand <<< "$row"
O=$W/out/$COH/$SRR; [ -f $O/DONE ] && { echo "$SRR already done"; exit 0; }
mkdir -p $O $W/tmp/$SRR; D=$W/tmp/$SRR; FQ=$D/R1.fastq.gz
mkdir -p /home/user/e3/results      # bench.py (from e3_cost/scripts) hard-codes this output folder
case $COH in munger) LAY=100;; *) LAY=50;; esac
if [ $strand = rf ]; then KS="--rf-stranded"; SL=SR; FS=2; else KS=""; SL=U; FS=0; fi
BENCH="$T/python $S/bench.py $LAY $COH/$SRR"; REP=1
b(){ local m=$1; shift; $BENCH $m $TH $REP -- "$@" | tail -1 >> $O/timing.tsv; }     # row: layout sample method threads rep wall cpu rss_mb rc
log(){ echo "$(date -u +%H:%M:%S) $SRR $*"; }
fail(){ log "FAILED $1"; echo "$1" > $O/FAILED; exit 1; }
rm -f $O/FAILED $O/timing.tsv
log "start $COH $title reads=$nreads strand=$strand"
# 1. download (not timed)
for i in 1 2 3; do
  $T/aria2c -q -x 8 -s 8 --allow-overwrite=true -d $D -o R1.fastq.gz "https://$ftp" && break; sleep $((i*20))
done
[ -s $FQ ] || fail download
ls -l $FQ | awk '{print $5}' > $O/fastq_bytes.txt
# 2. kallisto ref, pan
KPAN=${KPAN:-munger}; KCAP=${KPAN_CAP:-0}     # kallisto-pan runs on: all | munger (Munger only) | off ; KPAN_CAP = seconds cap per sample (0 = none)
for m in ref pan; do
  if [ $m = pan ]; then
    if [ $KPAN = off ] || { [ $KPAN = munger ] && [ $COH != munger ]; }; then echo "skipped by policy KPAN=$KPAN" > $O/kallisto_pan.SKIPPED; continue; fi
  fi
  CAPCMD=""; [ $m = pan ] && [ $KCAP -gt 0 ] && CAPCMD="timeout $KCAP "
  b kallisto_$m "${CAPCMD}$T/kallisto quant -i $W/idx/kallisto_$m.idx -o $D/kal_$m --single -l 250 -s 25 $KS -t $TH --plaintext $FQ > $D/kal_$m.log 2>&1"
  if [ ! -s $D/kal_$m/abundance.tsv ]; then
    if [ $m = pan ] && [ $KCAP -gt 0 ]; then echo "capped at ${KCAP}s" > $O/kallisto_pan.TIMEOUT; continue; fi
    fail kallisto_$m
  fi
  $T/python $S/fcsum.py kallisto $D/kal_$m $W/fa/tx2gene.tsv $O/kallisto_$m.genes.tsv.gz > /dev/null
  cp $D/kal_$m/run_info.json $O/kallisto_${m}_run_info.json
done
# 3. countTags (single-threaded) on the Kmerator tag list
b counttags "$CT -k 31 -i $W/km/tags/tags.u.raw.gz --summary $O/counttags_summary.txt $FQ > $D/ct.tsv 2> $D/ct.err"
[ -s $D/ct.tsv ] || fail counttags
$T/python $S/fcsum.py counttags $D/ct.tsv $W/km/tags $O/kmerator.genes.tsv.gz; rm -f $D/ct.tsv
# 4. STAR + featureCounts
b star "$T/STAR --runThreadN $TH --genomeDir $W/idx/star_sparse --readFilesIn $FQ --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $D/star_ > $D/star.out 2>&1"
[ -s $D/star_Aligned.out.bam ] || fail star
cp $D/star_Log.final.out $O/star_Log.final.out
b featurecounts "$T/featureCounts -T $TH -s $FS -t exon -g gene_id -a $W/ref/genes.gtf -o $D/fc.txt $D/star_Aligned.out.bam > $D/fc.log 2>&1"
[ -s $D/fc.txt ] || fail featurecounts
$T/python $S/fcsum.py fcounts $D/fc.txt $O/featurecounts.genes.tsv.gz; cp $D/fc.txt.summary $O/featurecounts.summary.txt
rm -f $D/star_Aligned.out.bam
# 5. reduced-read Salmon: seqtk, seeds 20261008 (25%) and 20261009 (10%), drawn directly from the full file
for pc in 25 10; do
  if [ $pc = 25 ]; then SEED=20261008; FR=0.25; else SEED=20261009; FR=0.10; fi
  b seqtk_p$pc "$T/seqtk sample -s $SEED $FQ $FR > $D/sub$pc.fastq"
  [ -s $D/sub$pc.fastq ] || fail seqtk$pc
  echo "$pc $(( $(wc -l < $D/sub$pc.fastq) / 4 )) $nreads" >> $O/subsample_reads.txt
  for m in ref pan; do
    b salmon_${m}_p$pc "$SAL quant -i $W/idx/salmon_$m -l $SL -r $D/sub$pc.fastq -p $TH -o $D/sal_${m}_p$pc -q > $D/sal_${m}_p$pc.log 2>&1"
    [ -s $D/sal_${m}_p$pc/quant.sf ] || fail salmon_${m}_p$pc
    $T/python $S/fcsum.py salmon $D/sal_${m}_p$pc/quant.sf $W/fa/tx2gene.tsv $O/salmon_${m}_p$pc.genes.tsv.gz > /dev/null
    cp $D/sal_${m}_p$pc/aux_info/meta_info.json $O/salmon_${m}_p${pc}_meta_info.json
  done
  rm -f $D/sub$pc.fastq
done
# 6. integrity numbers: reads processed by each tool vs the ENA read_count
$T/python - <<PY > $O/integrity.json
import json
o='$O'; out={'ena_read_count':$nreads}
k=json.load(open(o+'/kallisto_ref_run_info.json')); out['kallisto_ref_n_processed']=k['n_processed']; out['kallisto_ref_n_pseudoaligned']=k['n_pseudoaligned']
import os
if os.path.exists(o+'/kallisto_pan_run_info.json'):
    k=json.load(open(o+'/kallisto_pan_run_info.json')); out['kallisto_pan_n_processed']=k['n_processed']; out['kallisto_pan_n_pseudoaligned']=k['n_pseudoaligned']
else: out['kallisto_pan']='skipped' if os.path.exists(o+'/kallisto_pan.SKIPPED') else 'timeout'
for l in open(o+'/star_Log.final.out'):
    if 'Number of input reads' in l: out['star_input_reads']=int(l.split('|')[1])
    if 'Uniquely mapped reads number' in l: out['star_unique']=int(l.split('|')[1])
for l in open(o+'/counttags_summary.txt'): 
    if 'total_reads' in l: out['counttags_total_reads_line']=l.strip()
for s in ('ref','pan'):
  for pc in (10,25):
    m=json.load(open(f'{o}/salmon_{s}_p{pc}_meta_info.json')); out[f'salmon_{s}_p{pc}_processed']=m['num_processed']; out[f'salmon_{s}_p{pc}_mapped']=m['num_mapped']
print(json.dumps(out,indent=1))
PY
rm -rf $D; touch $O/DONE; log "done"
