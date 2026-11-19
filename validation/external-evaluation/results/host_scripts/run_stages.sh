#!/bin/bash
# Stages of validation/run_external_evaluation.py on the host, each in the conda
# environment that holds its tools.
# Usage: bash run_stages.sh "<stages>" "<isolates, comma-separated, or empty for all>" <threads> [extra runner arguments]
set -eo pipefail
STAGES="$1"; ISOLATES="$2"; THREADS="$3"; shift 3
BASE=/data/PROJECTS/locus_recon_revision/external
SRC=/data/PROJECTS/locus_recon_revision/repo        # read-only checkout at the protocol commit
REPO=$BASE/repo_run                                 # working copy that receives the runner
WORK=$BASE/work
OUT=$BASE/results/external-evaluation-results
JOBDIR=$(pwd)
source /home/crodsa/conda/etc/profile.d/conda.sh
mkdir -p "$BASE/results" "$WORK" "$BASE/logs"
if [ ! -d "$REPO" ]; then cp -r "$SRC" "$REPO"; fi
if ! cmp -s "$JOBDIR/run_external_evaluation.py" "$REPO/validation/run_external_evaluation.py"; then
  cp "$JOBDIR/run_external_evaluation.py" "$REPO/validation/run_external_evaluation.py"
fi
cd "$REPO"
echo "repo copy at $(git rev-parse HEAD); runner sha256 $(sha256sum validation/run_external_evaluation.py | cut -d' ' -f1)"
SEL=""; if [ -n "$ISOLATES" ]; then SEL="--isolates $ISOLATES"; fi
for STAGE in $STAGES; do
  if [ "$STAGE" = comparators ] || [ "$STAGE" = concordance ]; then conda activate lr-comparators; else conda activate locus-recon-101; fi
  echo "=== $STAGE $(date -u +%FT%TZ) env $CONDA_DEFAULT_ENV isolates ${ISOLATES:-all}"
  python validation/run_external_evaluation.py --workdir "$WORK" --outdir "$OUT" --stage "$STAGE" \
      --threads "$THREADS" --seqtk /home/crodsa/conda/envs/bioinfo/bin/seqtk $SEL "$@" 2>&1 \
      | tee -a "$BASE/logs/${STAGE}.log"
  conda deactivate
done
cd "$JOBDIR"
case " $STAGES " in *" score "*)
  tar czf external-evaluation-results.tar.gz -C "$BASE/results" external-evaluation-results
  cp "$OUT"/external_summary.tsv "$OUT"/decision_tables.tsv . 2>/dev/null || true
  ls -lh external-evaluation-results.tar.gz || true;;
esac
echo "=== done $(date -u +%FT%TZ)"
