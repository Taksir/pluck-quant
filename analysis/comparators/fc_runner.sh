#!/bin/bash
# Cloud runner K (1..6) of the accuracy comparison.  usage: bash fc_runner.sh K     (safe to re-run: resumes)
# Works in $W on local disk (ephemeral).  Slice = rows of scripts/runner_slices.tsv with runner == K.  After each sample: copy to final_comparison/out/COHORT/SRR/ and
# read the sha256 back from the shared folder (silent drops happen when the quota fills).  Progress -> runners/rK.md, resource rows -> timing/rK.tsv.
set -uo pipefail
# Folder constants below point to the original working folder; set them for your copy.
K=$1; F=/mnt/project-files/revision/final_comparison; SH=$F/scripts; export W=/home/user/fc; STOP_MB=8500
mkdir -p $F/runners $F/timing $F/out; PROG=$F/runners/r$K.md
log(){ echo "- $(date -u +%m-%dT%H:%MZ) $*" >> $PROG; echo "$*"; }
log "runner $K start"
mkdir -p $W; cd $W
if [ ! -f $W/BOOT_OK ]; then
  # the shared scripts may still be syncing: wait until every file in the manifest verifies
  for i in $(seq 1 40); do (cd $SH && sha256sum -c MANIFEST.sha256 --quiet >/dev/null 2>&1) && break; [ $i = 1 ] && log "runner $K: waiting for scripts to sync"; sleep 60; done
  (cd $SH && sha256sum -c MANIFEST.sha256 --quiet) || { log "runner $K: SCRIPT MANIFEST CHECK FAILED, stopping"; exit 1; }
  log "runner $K: scripts verified; boot (env, refs, founders, tags, indexes)"
  SH=$SH bash $SH/fc_boot.sh > $W/boot.out 2>&1 || { log "runner $K: BOOT FAILED: $(tail -3 $W/boot.out | tr '\n' ' ')"; exit 1; }
  log "runner $K: boot done"
fi
cp -n $SH/samples_fc.tsv $W/ 2>/dev/null; mkdir -p $W/scripts; cp -n $SH/*.sh $SH/*.py $W/scripts/ 2>/dev/null
# tag-file identity check against the shared tag hash
cmp -s $W/km/tags/tags.sha256.shared $SH/kmerator_tags/tags.sha256 || { log "runner $K: TAG HASH FILE DIFFERS from shared, stopping"; exit 1; }
used_ok(){ local u=$(du -sm /mnt/project-files 2>/dev/null | cut -f1); [ "${u:-0}" -lt $STOP_MB ] || { log "runner $K: shared folder at ${u} MB >= $STOP_MB stop line, stopping"; return 1; }; }
shared_ok(){ [ -s $F/out/$1/$2/SUMS.sha256 ] && (cd $F/out/$1/$2 && sha256sum -c SUMS.sha256 --quiet >/dev/null 2>&1) && [ -f $F/out/$1/$2/DONE ]; }
copy_verify(){ local c=$1 s=$2 try; local src=$W/out/$c/$s dst=$F/out/$c/$s
  for try in 1 2 3; do
    mkdir -p $dst; (cd $src && sha256sum $(ls | grep -v '^DONE$') > /tmp/sums.$s) || return 1
    cp $src/* $dst/ 2>/dev/null; cp /tmp/sums.$s $dst/SUMS.sha256; sleep 10
    (cd $dst && sha256sum -c SUMS.sha256 --quiet >/dev/null 2>&1) && { touch $dst/DONE; shared_ok $c $s && return 0; }
    log "runner $K: $c/$s copy not verified (try $try)"; sleep 60
  done; return 1; }
n=0; unverified=0; failed=0
while IFS=$'\t' read -r coh srr rk est; do
  [ "$rk" = "$K" ] || continue
  shared_ok $coh $srr && continue
  n=$((n+1)); [ $((n % 10)) = 1 ] && { used_ok || exit 3; }
  ok=0
  for try in 1 2; do bash $W/scripts/fc_sample.sh $coh $srr > $W/logs/sample_$srr.log 2>&1 && ok=1 && break; log "runner $K: $coh/$srr failed try $try ($(tail -1 $W/logs/sample_$srr.log))"; rm -rf $W/tmp/$srr; done
  [ $ok = 1 ] || { log "runner $K: $coh/$srr FAILED twice, skipped"; failed=$((failed+1)); continue; }
  cat $W/out/$coh/$srr/timing.tsv >> $F/timing/r$K.tsv 2>/dev/null
  if copy_verify $coh $srr; then log "runner $K: $coh/$srr done (copy verified)"; else unverified=$((unverified+1)); log "runner $K: $coh/$srr computed, COPY NOT VERIFIED (kept locally; re-run)"; fi
done < <(tail -n +2 $SH/runner_slices.tsv)
log "runner $K: finished; failed=$failed unverified=$unverified"
[ $failed = 0 ] && [ $unverified = 0 ] && log "runner $K: ALL DONE" || log "runner $K: INCOMPLETE (re-run the same command)"
