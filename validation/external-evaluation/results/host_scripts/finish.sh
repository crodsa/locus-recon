#!/bin/bash
# Last steps on the host: ARIBA and KMA for every isolate, a second subsampling to
# show that the read subsets are reproduced, scoring, and a second independent run
# of one isolate with the final runner.
set -eo pipefail
STEPS="${1:-comparators concordance verify score rerun pack}"
BASE=/data/PROJECTS/locus_recon_revision/external
OUT=$BASE/results/external-evaluation-results
RERUN=FDAARGOS_262
JOBDIR=$(pwd)
source /home/crodsa/conda/etc/profile.d/conda.sh
for STEP in $STEPS; do
  cd "$JOBDIR"
  case $STEP in
    comparators) bash ./run_stages.sh "comparators" "" 12;;
    verify) bash ./run_stages.sh "prepare" "" 12 --verify-subsample;;
    concordance) bash ./run_stages.sh "concordance" "" 12;;
    score) bash ./run_stages.sh "score" "" 12;;
    rerun)
      # Same commands as the main run, in an empty working directory; only the
      # downloaded files are shared.
      mkdir -p $BASE/work_rerun
      [ -e $BASE/work_rerun/downloads ] || ln -s $BASE/work/downloads $BASE/work_rerun/downloads
      cd $BASE/repo_run
      for STAGE in prepare assemble reconstruct copy-number block-size comparators concordance score; do
        if [ "$STAGE" = comparators ] || [ "$STAGE" = concordance ]; then conda activate lr-comparators; EXTRA="--threads-per-run 4"; else conda activate locus-recon-101; EXTRA="--assembly-jobs 1 --assembly-threads 8"; fi
        python validation/run_external_evaluation.py --workdir $BASE/work_rerun --outdir $BASE/results/rerun_check \
          --stage $STAGE --isolates $RERUN --threads 8 --seqtk /home/crodsa/conda/envs/bioinfo/bin/seqtk $EXTRA 2>&1 | tee -a $BASE/logs/rerun_$STAGE.log | tail -n 3
        conda deactivate
      done
      conda activate locus-recon-101
      python "$JOBDIR/compare_rerun.py" $OUT $BASE/results/rerun_check $RERUN $OUT/rerun_check_$RERUN.tsv
      conda deactivate;;
    pack)
      cd "$JOBDIR"
      for NAME in README.md DEVIATIONS.md; do if [ -f "$JOBDIR/$NAME" ]; then cp "$JOBDIR/$NAME" $OUT/$NAME; fi; done
      mkdir -p $OUT/host_scripts $OUT/logs $OUT/concordance_reports
      for NAME in run_stages.sh follow.sh finish.sh compare_rerun.py; do if [ -f "$JOBDIR/$NAME" ]; then cp "$JOBDIR/$NAME" $OUT/host_scripts/; fi; done
      for NAME in host_probe_2026-11-18.txt bait_symbol_probe_2026-11-18.txt; do if [ -f "$JOBDIR/$NAME" ]; then cp "$JOBDIR/$NAME" $OUT/logs/; fi; done
      cp $BASE/logs/*.log $OUT/logs/ && gzip -n -f $OUT/logs/*.log
      cp $BASE/work/concordance/*.report $OUT/concordance_reports/
      for UNIT in $OUT/per_unit/*/*; do
        cp $BASE/work/isolates/$(basename $(dirname $UNIT))/truth/$(basename $UNIT).blastn.tsv $UNIT/bait_vs_complete_genome.blastn.tsv
      done
      (cd $OUT && find . -type f ! -name SHA256SUMS | sort | xargs sha256sum > SHA256SUMS)
      tar czf external-evaluation-results.tar.gz -C $BASE/results external-evaluation-results
      mkdir -p tables; cp $OUT/*.tsv $OUT/PROVENANCE.json tables/
      tar tzf external-evaluation-results.tar.gz | wc -l
      ls -lh external-evaluation-results.tar.gz
      du -sh $BASE/work $BASE/work_rerun $BASE/results 2>/dev/null || true;;
  esac
done
echo "=== finish.sh done $(date -u +%FT%TZ)"
