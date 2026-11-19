#!/bin/bash
# Process each isolate as soon as the prepare stage has written its reads.
# Usage: bash follow.sh "<isolates in order, space-separated>" <threads>
set -eo pipefail
ISOLATES="$1"; THREADS="$2"; STAGES="${3:-assemble reconstruct copy-number block-size comparators}"
WORK=/data/PROJECTS/locus_recon_revision/external/work
# One follower at a time works on an isolate and a group of stages: the claim is a
# directory, which only one process can create.
KIND=main; if [ "$STAGES" = comparators ]; then KIND=comparators; fi
mkdir -p "$WORK/claims"
for ISOLATE in $ISOLATES; do
  RECORD=$WORK/prepare/isolates/$ISOLATE.json
  until [ -f "$RECORD" ] && grep -q '"pairs_kept"' "$RECORD"; do sleep 60; done
  if ! mkdir "$WORK/claims/$KIND.$ISOLATE" 2>/dev/null; then
    echo "##### $ISOLATE ($KIND) is claimed by another follower; skipped"; continue
  fi
  echo "##### $ISOLATE ($KIND) $(date -u +%FT%TZ)"
  for STAGE in $STAGES; do
    if [ "$STAGE" = comparators ]; then
      bash ./run_stages.sh "comparators" "$ISOLATE" "$THREADS" --threads-per-run "$THREADS"
    else
      bash ./run_stages.sh "$STAGE" "$ISOLATE" "$THREADS" --assembly-jobs 1 --assembly-threads 8
    fi
  done
  touch "$WORK/claims/$KIND.$ISOLATE/done"
done
echo "##### all isolates done $(date -u +%FT%TZ)"
