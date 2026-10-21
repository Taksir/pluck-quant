#!/bin/bash
# Runner/check bootstrap on a cloud VM (4 vCPU, 15 GB): environment, references, founder transcriptomes, indexes, Kmerator tag files.  Idempotent; ends with $W/BOOT_OK.
# usage: SH=<shared scripts dir> bash fc_boot.sh      (W defaults to the original working folder; set it for your copy)
set -uo pipefail
# Folder constants below point to the original working folder; set them for your copy.
W=${W:-/home/user/fc}; SH=${SH:-/mnt/project-files/revision/final_comparison/scripts}; export W
mkdir -p $W/{scripts,code,logs,ref/vcf,fa,idx,km,tmp,out}; cd $W
[ -f BOOT_OK ] && { echo "boot already done"; exit 0; }
cp -n $SH/*.sh $SH/*.py $W/scripts/ 2>/dev/null; cp -n $SH/e9.py $SH/build_founders.py $W/code/ 2>/dev/null; cp -n $SH/samples_fc.tsv $W/
log(){ echo "$(date -u +%H:%M:%S) $*"; }
T=$W/mm/envs/tools/bin
log "env"; bash $SH/s1_env.sh > logs/s1_env.log 2>&1 || { log "ENV FAILED"; exit 1; }
if [ ! -x $W/src/countTags/bin/countTags ]; then
  mkdir -p $W/src; rm -rf $W/src/countTags; git clone -q https://github.com/Transipedia/countTags $W/src/countTags && (cd $W/src/countTags && git checkout -q 837b8e49258646ea77ca3644b62c5d407041584e && make >/dev/null 2>&1) || { log "countTags FAILED"; exit 1; }
fi
log "refs+founders"; bash $SH/s2_refs.sh > logs/s2_refs.log 2>&1; [ -s fa/stats.json ] || { log "FOUNDERS FAILED"; tail -5 logs/s2_refs.log; exit 1; }
rm -rf ref/vcf fa/vars_*.json
log "tags"; mkdir -p km/tags
$T/python scripts/mktags.py make $SH/kmerator_tags/contigs.fa.gz $SH/kmerator_tags km/tags > logs/mktags_make.log 2>&1 || { log "TAGS FAILED (hash check?)"; tail -3 logs/mktags_make.log; exit 1; }
cp $SH/kmerator_tags/tags.sha256 km/tags/tags.sha256.shared
$T/python scripts/mkdummy.py km/tags/dummy.fq && $W/src/countTags/bin/countTags -k 31 -i km/tags/tags.u.raw.gz --summary km/tags/dummy.sum km/tags/dummy.fq > km/tags/dummy.tsv 2>/dev/null \
  && $T/python scripts/fcsum.py ctperm km/tags/dummy.tsv km/tags/tags.u.raw.gz km/tags > logs/ctperm.log 2>&1 && rm -f km/tags/dummy.tsv || { log "CTPERM FAILED"; exit 1; }
cmp -s km/tags/ct_order.sha256 $SH/kmerator_tags/ct_order.sha256 || { log "CT ORDER DIFFERS FROM REFERENCE (countTags output order)"; exit 1; }
log "indexes"; bash $SH/s3_idx.sh > logs/s3.log 2>&1; tail -1 logs/s3.log | grep -q S3_IDX_DONE || { log "INDEX FAILED"; tail -5 logs/s3.log; exit 1; }
rm -f ref/dna.fa ref/dna.fa.gz; touch BOOT_OK; log "BOOT_OK"; df -h $W | tail -1
