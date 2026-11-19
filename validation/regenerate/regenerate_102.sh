#!/usr/bin/env bash
# Rerun with Locus-Recon 1.0.2 the validation analyses that release 1.0.2 can
# change, and compare them with the records of release 1.0.1.
#
# Release 1.0.2 differs from 1.0.1 in locus_recon/copy_number.py (one rule of
# reconcile_copy_number) and in the version string.  Only these analyses are
# rerun: the 23S rRNA gene and gyrB of the five H. pylori isolates (its runner
# calls reconcile_copy_number), the standard workflow of LIBA-6656 with its
# release verification, locus-recon-copy-number without a graph on the
# LIBA-6656 tcdB alignment, and the simulated benchmark of 13 samples.
#
# The working directory of the 1.0.1 regeneration ($PREV) is read, never
# written: reads, the draft of LIBA-6656, the five SPAdes drafts of the
# H. pylori isolates and the LIBA-6656 alignment are taken from it by copy or
# symbolic link.  Each step writes $REGEN/logs/<step>.log (with the commands),
# a line in $REGEN/logs/steps.tsv and a marker in $REGEN/done/.
#
# Usage:
#   bash regenerate_102.sh prep            # clone at the commit, check digest and versions, inputs
#   bash regenerate_102.sh chain_run       # the four analyses (at most 6 threads)
#   bash regenerate_102.sh chain_package   # evidence, records, tests, comparison, package
#   bash regenerate_102.sh run <step function> [arguments]

REGEN=${REGEN:-/data/PROJECTS/locus_recon_revision/regen102}
PREV=${PREV:-/data/PROJECTS/locus_recon_revision/regen}
GIT_URL=${GIT_URL:-https://github.com/crodsa/locus-recon}
EXPECTED_COMMIT=${EXPECTED_COMMIT:-f985906673912a900fd0035e3045cd2a5f9a29f3}
PREVIOUS_COMMIT=${PREVIOUS_COMMIT:-3669742a69f57701c742d5a98f1037626df39c12}
EXPECTED_DIGEST=${EXPECTED_DIGEST:-a24d41d80d2a3e04729eda04efa6938d82ac2d99e4ea875377353f014aae9b35}
EXPECTED_VERSION=${EXPECTED_VERSION:-1.0.2}
PREVIOUS_VERSION=${PREVIOUS_VERSION:-1.0.1}
READ_CACHE=${READ_CACHE:-/home/crodsa/claude_scratch/locus_recon_comparators/downloads}
CONDA_SH=${CONDA_SH:-/home/crodsa/conda/etc/profile.d/conda.sh}
CONDA_ENV=${CONDA_ENV:-locus-recon-101}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$REGEN/repo
WORK=$REGEN/work
LOGS=$REGEN/logs
DONE=$REGEN/done
INPUTS=$REGEN/inputs
RECORDS=$REGEN/records
SCRIPTS=$REGEN/scripts
PACKAGE=$REGEN/package/regenerated-$EXPECTED_VERSION
PYDEPS=$REGEN/pydeps
HP_RUNS="SRR18335541 SRR18335539 SRR18335523 SRR18335501 SRR18335507"
STRAINS="Hpfe0001 Hpfe0002 Hpfe0003 Hpfe0004 Hpfe0006"

case "$REGEN" in
    "$PREV"|"$PREV"/*) echo "REGEN must not be inside the directory of the 1.0.1 run"; exit 2 ;;
esac

# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$CONDA_ENV"
set -uo pipefail
mkdir -p "$WORK" "$LOGS" "$DONE" "$INPUTS" "$RECORDS" "$SCRIPTS"
# The environment holds an installed copy of release 1.0.1.  The clone of
# 1.0.2 is put first on the module path; prep stops unless every command then
# reports 1.0.2, and every output records the version that wrote it.
export PYTHONPATH=$REPO
export PYTHONDONTWRITEBYTECODE=1

step() {
    local name=$1; shift
    if [ -e "$DONE/$name" ]; then
        echo "[skip] $name ($(cat "$DONE/$name"))"
        return 0
    fi
    echo "[start] $name $(date -u +%FT%TZ)"
    local started; started=$(date +%s)
    ( set -exo pipefail; "$@" ) > "$LOGS/$name.log" 2>&1
    local rc=$?
    local line; line="$name	rc=$rc	seconds=$(( $(date +%s) - started ))	finished=$(date -u +%FT%TZ)"
    echo "$line" >> "$LOGS/steps.tsv"
    echo "$line" > "$DONE/$name"
    echo "[end] $line"
    return 0
}

link_reads() {  # link_reads <directory> <run>...
    local directory=$1; shift
    mkdir -p "$directory"
    local run mate
    for run in "$@"; do
        for mate in 1 2; do
            if [ -s "$READ_CACHE/${run}_${mate}.fastq.gz" ] && [ ! -e "$directory/${run}_${mate}.fastq.gz" ]; then
                ln -s "$READ_CACHE/${run}_${mate}.fastq.gz" "$directory/${run}_${mate}.fastq.gz"
            fi
        done
    done
}

sha_at() {  # sha_at <commit> <path>: SHA-256 of a file at a commit, or "absent"
    if git cat-file -e "$1:$2" 2>/dev/null; then
        git show "$1:$2" | sha256sum | cut -d' ' -f1
    else
        echo absent
    fi
}

# ----------------------------------------------------------------------- prep
do_prep() {
    [ -d "$REPO/.git" ] || git clone -q "$GIT_URL" "$REPO"
    cd "$REPO"
    git checkout -q --detach "$EXPECTED_COMMIT"
    git rev-parse HEAD
    [ "$(git rev-parse HEAD)" = "$EXPECTED_COMMIT" ] || { echo "UNEXPECTED COMMIT"; exit 3; }
    [ -z "$(git status --porcelain)" ] || { echo "TREE NOT CLEAN"; git status --porcelain; exit 3; }
    local digest
    digest=$(python -c "from validation.common import source_sha256; print(source_sha256('.'))")
    echo "package digest of the clone: $digest"
    [ "$digest" = "$EXPECTED_DIGEST" ] || { echo "DIGEST MISMATCH"; exit 3; }

    # which code runs, and what the four commands report
    cd "$REGEN"
    python -c "import locus_recon; print('module in use:', locus_recon.__file__, locus_recon.VERSION)"
    ( cd / && env -u PYTHONPATH python -c "import locus_recon; print('module installed in the environment, not used:', locus_recon.__file__, locus_recon.VERSION)" ) || true
    : > "$REGEN/COMMAND_VERSIONS.tsv"
    local command reported
    for command in locus-recon locus-recon-depth-ratio locus-recon-copy-number locus-recon-graph-paths; do
        reported=$("$command" --version 2>&1 | head -1)
        printf '%s\t%s\n' "$command" "$reported" >> "$REGEN/COMMAND_VERSIONS.tsv"
        case "$reported" in
            *"$EXPECTED_VERSION"*) ;;
            *) echo "UNEXPECTED VERSION: $command reports $reported"; exit 3 ;;
        esac
    done
    cat "$REGEN/COMMAND_VERSIONS.tsv"
    {
        echo "locus-recon	$(locus-recon --version 2>&1 | head -1)"
        echo "python	$(python --version 2>&1)"
        echo "spades	$(spades.py --version 2>&1 | head -1)"
        echo "bwa-mem2	$(bwa-mem2 version 2>&1 | tail -1)"
        echo "samtools	$(samtools --version | head -1)"
        echo "blastn	$(blastn -version | head -1)"
        echo "makeblastdb	$(makeblastdb -version | head -1)"
        echo "biopython	$(python -c 'import Bio; print(Bio.__version__)')"
        echo "numpy	$(python -c 'import numpy; print(numpy.__version__)')"
        echo "git	$(git --version)"
        echo "platform	$(uname -srm)"
        echo "cpu	$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
    } > "$REGEN/TOOL_VERSIONS.tsv"
    cat "$REGEN/TOOL_VERSIONS.tsv"
    conda list --explicit > "$REGEN/conda_explicit_${CONDA_ENV}.txt" || true

    # what differs between the two releases
    cd "$REPO"
    local table=$REGEN/package_files_${PREVIOUS_COMMIT:0:7}_vs_${EXPECTED_COMMIT:0:7}.tsv
    printf 'file\tsha256_at_%s\tsha256_at_%s\tsame\n' "${PREVIOUS_COMMIT:0:7}" "${EXPECTED_COMMIT:0:7}" > "$table"
    local file before after
    for file in $( { git ls-tree -r --name-only "$PREVIOUS_COMMIT" locus_recon
                     git ls-tree -r --name-only "$EXPECTED_COMMIT" locus_recon; } | sort -u ); do
        before=$(sha_at "$PREVIOUS_COMMIT" "$file")
        after=$(sha_at "$EXPECTED_COMMIT" "$file")
        printf '%s\t%s\t%s\t%s\n' "$file" "$before" "$after" \
            "$([ "$before" = "$after" ] && echo yes || echo no)" >> "$table"
    done
    cat "$table"
    git diff --stat "$PREVIOUS_COMMIT" "$EXPECTED_COMMIT" > "$REGEN/git_diff_stat_${PREVIOUS_COMMIT:0:7}_${EXPECTED_COMMIT:0:7}.txt"
    cat "$REGEN/git_diff_stat_${PREVIOUS_COMMIT:0:7}_${EXPECTED_COMMIT:0:7}.txt"
    git diff "$PREVIOUS_COMMIT" "$EXPECTED_COMMIT" -- locus_recon > "$REGEN/locus_recon_${PREVIOUS_VERSION}_to_${EXPECTED_VERSION}.diff"
    # every place of the package and of the validation scripts that reaches
    # reconcile_copy_number, directly or through the command that calls it
    grep -rnE "reconcile_copy_number|locus-recon-copy-number|copy_number_cli|locus_recon\.copy_number\b|from \.copy_number\b" \
        --include='*.py' locus_recon validation > "$REGEN/reconcile_copy_number_references.txt" || true
    cut -d: -f1 "$REGEN/reconcile_copy_number_references.txt" | sort | uniq -c
    ls validation
    git tag

    # inputs, from the directory of the 1.0.1 run
    mkdir -p "$INPUTS/reads"
    [ -s "$INPUTS/LIBA6656_ST154.fasta" ] || cp "$PREV/inputs/LIBA6656_ST154.fasta" "$INPUTS/LIBA6656_ST154.fasta"
    link_reads "$INPUTS/reads" $HP_RUNS ERR467623
    ls -l "$INPUTS/reads"
    [ -d "$PYDEPS" ] || cp -r "$PREV/pydeps" "$PYDEPS"
    env PYTHONPATH="$PYDEPS/pytest" python -c "import pytest; print('pytest', pytest.__version__)"
}

do_input_checksums() {
    cd "$INPUTS"
    : > INPUT_SHA256.tsv
    local path
    for path in reads/*.fastq.gz LIBA6656_ST154.fasta; do
        printf '%s\t%s\t%s\n' "$path" "$(stat -L -c %s "$path")" "$(sha256sum "$path" | cut -d' ' -f1)" >> INPUT_SHA256.tsv
    done
    cat INPUT_SHA256.tsv
    if diff <(sort INPUT_SHA256.tsv) <(sort "$PREV/inputs/INPUT_SHA256.tsv"); then
        echo "INPUTS_SAME_AS_THE_1.0.1_RUN yes"
    else
        echo "INPUTS_SAME_AS_THE_1.0.1_RUN no"
    fi
}

# ----------------------------------------------------------------------- mock
do_mock() {  # do_mock <results directory name>
    cd "$REPO"
    python validation/run_mock_validation.py \
        --results-dir "$WORK/mock/$1" --threads 2 --memory-gb 4 --force
}

do_reference_results() {
    cd "$REPO"
    python "$SCRIPTS/make_reference_results.py" \
        --results-dir "$WORK/mock/mock_results" --out "$RECORDS/reference_results"
    ls -l "$RECORDS/reference_results"
}

# --------------------------------------------------------------- low-copy 23S
do_low_copy() {
    cd "$REPO"
    link_reads "$WORK/low_copy_23S/reads" $HP_RUNS
    # The five SPAdes --isolate drafts are those of the 1.0.1 run (made by the
    # closed-genome stage of its tier calibration and reused by its 23S
    # analysis): the runner keeps a draft that already exists.  Reads are
    # mapped again and every later stage is run again.
    # The GenBank records of the closed genomes and of strain 26695 are those
    # the 1.0.1 run downloaded: the runner keeps a record that already exists,
    # so both runs audit the same truth files.
    local strain name
    for name in closed_truth fixed_baits; do
        rm -rf "$WORK/low_copy_23S/$name"
        cp -r "$PREV/work/low_copy_23S/$name" "$WORK/low_copy_23S/$name"
    done
    ( cd "$WORK/low_copy_23S" && sha256sum closed_truth/* fixed_baits/* )
    rm -rf "$WORK/low_copy_23S/results"
    : > "$WORK/low_copy_23S/drafts_reused.tsv"
    for strain in $STRAINS; do
        rm -rf "$WORK/low_copy_23S/samples/$strain"
        mkdir -p "$WORK/low_copy_23S/samples/$strain/spades"
        for name in scaffolds.fasta contigs.fasta assembly_graph_after_simplification.gfa params.txt spades.log; do
            cp "$PREV/work/tier/$strain/spades/$name" "$WORK/low_copy_23S/samples/$strain/spades/$name"
        done
        for name in contigs.fasta assembly_graph_after_simplification.gfa; do
            printf '%s\t%s\t%s\t%s\n' "$strain" "$name" \
                "$(sha256sum "$WORK/low_copy_23S/samples/$strain/spades/$name" | cut -d' ' -f1)" \
                "$(sha256sum "$PREV/work/low_copy_23S/samples/$strain/spades/$name" | cut -d' ' -f1)" \
                >> "$WORK/low_copy_23S/drafts_reused.tsv"
        done
    done
    cat "$WORK/low_copy_23S/drafts_reused.tsv"
    set +e
    python validation/workflows/run_low_copy_23S.py \
        --work-dir "$WORK/low_copy_23S" --threads 4 --memory-gb 8
    local rc=$?
    set -e
    # exit 1 is the outcome when a primary criterion is not met; evaluation.json
    # names the criterion
    echo "RUN_LOW_COPY_23S_EXIT $rc"
}

# ------------------------------------------------- LIBA-6656 standard workflow
do_liba_gfa() {  # do_liba_gfa <threads>
    local threads=$1
    local run=$WORK/liba6656_gfa_rerun/run_${threads}_threads
    rm -rf "$run"
    mkdir -p "$run"
    cd "$run"
    cp "$INPUTS/LIBA6656_ST154.fasta" LIBA6656_ST154.fasta
    printf 'sample_id\tassembly_path\tr1_path\tr2_path\n' > samples.tsv
    printf 'LIBA6656_ST154\t%s\t%s\t%s\n' "$run/LIBA6656_ST154.fasta" \
        "$INPUTS/reads/ERR467623_1.fastq.gz" "$INPUTS/reads/ERR467623_2.fastq.gz" >> samples.tsv
    python -m locus_recon.cli \
        --samplesheet samples.tsv \
        --main-output-dir output \
        --bait "$REPO/validation/liba6656-gfa-rerun/tcdB_111_ungapped_baits.fasta" \
        --locus tcdB \
        --threads "$threads" \
        --parallel-samples 1 \
        --memory-per-sample 12 \
        --min-discovery-identity 80 \
        --min-hsp-length 50 \
        --region-padding 300 \
        --no-progress
}

do_liba_verify() {  # do_liba_verify <threads of first run> <threads of second run>
    local base=$WORK/liba6656_gfa_rerun
    cd "$base"
    python -m locus_recon.graph_cli \
        --gfa "run_$1_threads/output/LIBA6656_ST154/spades_local/assembly_graph_after_simplification.gfa" \
        --output tcdB_terminal_paths_8000_9000.fasta \
        --summary tcdB_terminal_paths_8000_9000.tsv \
        --min-length 8000 \
        --max-length 9000 \
        --prefix LIBA6656_tcdB_terminal_path
    cd "$REPO"
    python validation/liba6656-gfa-rerun/verify_release_rerun.py \
        --run-dir "$base/run_$1_threads/output" --run-dir "$base/run_$2_threads/output" \
        --terminal-paths "$base/tcdB_terminal_paths_8000_9000.fasta" \
        --output "$base/RELEASE_VERIFICATION.json"
    if cmp "$base/tcdB_terminal_paths_8000_9000.fasta" \
        validation/liba6656-gfa-rerun/tcdB_terminal_paths_8000_9000.fasta; then
        echo "TERMINAL_PATHS_FASTA byte-identical to the deposit"
    else
        echo "TERMINAL_PATHS_FASTA differs in bytes from the deposit"
    fi
    local name
    for name in tcdB_terminal_paths_8000_9000.fasta tcdB_terminal_paths_8000_9000.tsv; do
        if cmp -s "$base/$name" "$PREV/work/liba6656_gfa_rerun/$name"; then
            echo "$name byte-identical to the file of the 1.0.1 run"
        else
            echo "$name differs in bytes from the file of the 1.0.1 run"
        fi
    done
    [ -z "$(git status --porcelain)" ] && echo "WORKING_CLONE_CLEAN yes" || echo "WORKING_CLONE_CLEAN no"
}

# ------------------------- locus-recon-copy-number without a graph, LIBA-6656
do_consensus_nograph() {
    local out=$WORK/liba6656_consensus_nograph
    rm -rf "$out"
    mkdir -p "$out"
    cd "$out"
    # the alignment of ERR467623 to the draft of LIBA-6656 made by the 1.0.1 run
    ln -s "$PREV/work/liba6656_depth/results/ERR467623_vs_LIBA6656_draft.sorted.bam" .
    ln -s "$PREV/work/liba6656_depth/results/ERR467623_vs_LIBA6656_draft.sorted.bam.bai" .
    local arguments=() region start end compartment
    while IFS=$'\t' read -r region start end compartment; do
        [ "$region" = "assembly_region" ] && continue
        [ -z "$region" ] && continue
        arguments+=(--locus "$region" --gene-span "${start}-${end}")
    done < "$REPO/validation/manifests/liba6656_depth_regions.tsv"
    echo "regions: $(( ${#arguments[@]} / 4 ))"
    printf '%s ' locus-recon-copy-number --bam ERR467623_vs_LIBA6656_draft.sorted.bam "${arguments[@]}" \
        --gene-length 7104 --output-json consensus_without_graph.json \
        --output-tsv consensus_without_graph.tsv > command.txt
    echo >> command.txt
    set +e
    locus-recon-copy-number --bam "$out/ERR467623_vs_LIBA6656_draft.sorted.bam" "${arguments[@]}" \
        --gene-length 7104 --output-json consensus_without_graph.json \
        --output-tsv consensus_without_graph.tsv > stdout.txt 2> stderr.txt
    local rc=$?
    set -e
    echo "COPY_NUMBER_EXIT $rc"
    cat stdout.txt
    cat stderr.txt
    sha256sum "$PREV/work/liba6656_depth/results/ERR467623_vs_LIBA6656_draft.sorted.bam"
}

# -------------------------------------------------------------------- package
collect() {  # collect <command>
    python "$SCRIPTS/collect_102.py" "$1" --regen "$REGEN" --package "$PACKAGE" --repo "$REPO" \
        --previous "$PREV" --version "$EXPECTED_VERSION" --previous-version "$PREVIOUS_VERSION" \
        --digest "$EXPECTED_DIGEST" --commit "$EXPECTED_COMMIT" --previous-commit "$PREVIOUS_COMMIT" \
        --notes "$SCRIPTS/difference_notes_102.tsv"
}

# Static check of which runners can reach reconcile_copy_number (see the script).
do_reach() {
    python "$SCRIPTS/import_reach_102.py" "$REPO" "$REGEN/reconcile_copy_number_reachability.tsv"
    # every import line of the package and of the validation scripts that names copy_number
    grep -rnE "^\s*(from|import)\b.*copy_number" --include='*.py' "$REPO/locus_recon" "$REPO/validation" \
        | sed "s#$REPO/##" || true
    # Records of analyses that are not part of this rerun and that hold the
    # method to which the changed rule applies: per file, the lines that name
    # the method and, of these, the lines that also hold SUPPORTED.
    cd "$REPO"
    local table=$REGEN/records_not_rerun_that_hold_the_method.tsv file
    printf 'file\tlines_with_GRAPH_COUNT_OVER_DISCORDANT_DEPTH\tof_these_lines_with_SUPPORTED\n' > "$table"
    for file in $(grep -rl "GRAPH_COUNT_OVER_DISCORDANT_DEPTH" validation \
                      --include='*.tsv' --include='*.json' --include='*.csv' | sort); do
        case "$file" in validation/low-copy-23S/*) continue ;; esac
        printf '%s\t%s\t%s\n' "$file" "$(grep -c GRAPH_COUNT_OVER_DISCORDANT_DEPTH "$file")" \
            "$(grep GRAPH_COUNT_OVER_DISCORDANT_DEPTH "$file" | grep -c SUPPORTED || true)" >> "$table"
    done
    cat "$table"
    # The comparison reads the 1.0.1 records from validation/ of the tree.
    # Check that, for the folders rerun here, the tree holds the files of the
    # package of the 1.0.1 regeneration with the same bytes.
    local previous=$PREV/package/regenerated-$PREVIOUS_VERSION/records folder
    {
        echo "records of the $PREVIOUS_VERSION regeneration package against validation/ of the tree at $EXPECTED_COMMIT"
        for folder in low-copy-23S liba6656-gfa-rerun liba6656-depth reference_results; do
            echo "== $folder"
            diff -rq "$previous/$folder" "validation/$folder" | sed "s#$previous/##; s#$PREV/##" || true
            echo "files of the package that differ in the tree: $(diff -rq "$previous/$folder" "validation/$folder" | grep -c '^Files ' || true)"
        done
        echo "== INTEGRATION_RUN.md"
        cmp -s "$previous/INTEGRATION_RUN.md" validation/INTEGRATION_RUN.md && echo "same bytes" || echo "DIFFERENT"
    } > "$REGEN/baseline_records_check.txt"
    cat "$REGEN/baseline_records_check.txt"
}

do_evidence() { collect evidence; }
do_records() { collect records; }
do_integration() { collect integration; }
do_compare() { collect compare; }
do_finish() { collect finish; }
do_summary() { collect summary; }

# The test suite of the commit, as committed and with the regenerated records
# copied over validation/; a third tree applies the suggested test update when
# one is supplied.
do_tests() {
    local base=$WORK/tests tree
    rm -rf "$base" "$PACKAGE/tests"
    mkdir -p "$base" "$PACKAGE/tests"
    for tree in release_as_committed release_with_regenerated_records \
                release_with_regenerated_records_and_test_update; do
        if [ "$tree" = release_with_regenerated_records_and_test_update ] && \
           [ ! -s "$SCRIPTS/make_test_update_102.py" ]; then
            continue
        fi
        git clone -q "$REPO" "$base/$tree"
        ( cd "$base/$tree"
          git checkout -q --detach "$EXPECTED_COMMIT"
          case "$tree" in
              release_with_regenerated_records)
                  cp -r "$PACKAGE/records/." validation/ ;;
              release_with_regenerated_records_and_test_update)
                  # the update is written by a script, saved as a patch, undone,
                  # and applied again from the patch file: the file that is
                  # delivered is the one that is tested
                  python "$SCRIPTS/make_test_update_102.py"
                  git diff -- tests > "$base/suggested_test_update_1.0.2.patch"
                  git checkout -q -- tests
                  [ -z "$(git status --porcelain)" ]
                  git apply --verbose "$base/suggested_test_update_1.0.2.patch"
                  git diff --stat -- tests
                  cp -r "$PACKAGE/records/." validation/ ;;
          esac
          git rev-parse HEAD
          git status --porcelain
          env -u PYTHONPATH python -c "import locus_recon; print('package under test:', locus_recon.__file__, locus_recon.VERSION)"
          set +e
          env PYTHONPATH="$PYDEPS/pytest" python -m pytest -q -p no:cacheprovider \
              > "$base/$tree.pytest.txt" 2>&1
          echo "PYTEST_EXIT $tree $?"
          set -e
          tail -n 40 "$base/$tree.pytest.txt" )
        cp "$base/$tree.pytest.txt" "$PACKAGE/tests/"
    done
    if [ -s "$base/suggested_test_update_1.0.2.patch" ]; then
        cp "$base/suggested_test_update_1.0.2.patch" "$PACKAGE/tests/"
        cat "$base/suggested_test_update_1.0.2.patch"
    fi
}

do_manifest() {
    mkdir -p "$PACKAGE/scripts/validation/regenerate" "$PACKAGE/logs" "$PACKAGE/inputs"
    local name
    for name in regenerate_102.sh collect_102.py collect_regenerated.py make_reference_results.py \
                difference_notes_102.tsv make_test_update_102.py import_reach_102.py; do
        [ -s "$SCRIPTS/$name" ] && cp "$SCRIPTS/$name" "$PACKAGE/scripts/validation/regenerate/"
    done
    find "$LOGS" -maxdepth 1 -type f -size -5000k -exec cp {} "$PACKAGE/logs/" \;
    cp "$INPUTS/INPUT_SHA256.tsv" "$REGEN/TOOL_VERSIONS.tsv" "$REGEN/COMMAND_VERSIONS.tsv" \
       "$REGEN"/conda_explicit_*.txt "$REGEN"/package_files_*.tsv "$REGEN"/git_diff_stat_*.txt \
       "$REGEN"/locus_recon_*.diff "$REGEN/reconcile_copy_number_references.txt" \
       "$REGEN/reconcile_copy_number_reachability.tsv" \
       "$REGEN/records_not_rerun_that_hold_the_method.tsv" \
       "$REGEN/baseline_records_check.txt" "$PACKAGE/inputs/"
    for name in README.md DEVIATIONS.md; do
        if [ -s "$SCRIPTS/$name" ]; then cp "$SCRIPTS/$name" "$PACKAGE/$name"; else rm -f "$PACKAGE/$name"; fi
    done
    # the two tables of the README are written from the delivered files
    if [ -s "$PACKAGE/README.md" ]; then
        python - "$PACKAGE/README.md" "$REGEN"/package_files_*.tsv \
            "$REGEN/reconcile_copy_number_reachability.tsv" <<'EOF'
import csv
import sys

readme, files, reach = sys.argv[1:4]
text = open(readme).read()
rows = list(csv.reader(open(files), delimiter="\t"))
table = [f"| File | SHA-256 at `{rows[0][1].rsplit('_', 1)[-1]}` | SHA-256 at "
         f"`{rows[0][2].rsplit('_', 1)[-1]}` | Same |", "|---|---|---|---|"]
table += [f"| `{row[0]}` | `{row[1]}` | `{row[2]}` | {row[3]} |" for row in rows[1:]]
assert text.count("@@PACKAGE_FILES_TABLE@@") == 1
text = text.replace("@@PACKAGE_FILES_TABLE@@", "\n".join(table))
rows = list(csv.DictReader(open(reach), delimiter="\t"))
table = ["| Analysis | Rerun with 1.0.2 | Reaches the function by imports and commands | "
         "Modules on that path that refer to it | Only if a script named in text is followed |",
         "|---|---|---|---|---|"]
for row in rows:
    table.append(
        f"| {row['analysis']} | {row['rerun_with_1.0.2']} | "
        f"{row['can_reach_reconcile_copy_number_by_imports_and_commands']} | "
        f"{row['modules_reached_by_imports_and_commands_that_refer_to_reconcile_copy_number']} | "
        f"{row['further_referring_modules_if_scripts_named_in_text_are_followed']} |")
assert text.count("@@REACHABILITY_TABLE@@") == 1
text = text.replace("@@REACHABILITY_TABLE@@", "\n".join(table))
open(readme, "w").write(text)
EOF
    fi
    collect manifest
}

do_pack() {
    cd "$REGEN/package"
    tar -czf "$HERE/regenerated-$EXPECTED_VERSION.tar.gz" "regenerated-$EXPECTED_VERSION"
    local name
    for name in old_vs_new_${PREVIOUS_VERSION}_vs_${EXPECTED_VERSION}.tsv \
                old_vs_new_${PREVIOUS_VERSION}_vs_${EXPECTED_VERSION}_summary.tsv \
                summary_${EXPECTED_VERSION}.tsv PROVENANCE.json SHA256SUMS.txt FILES_LEFT_ON_HOST.tsv \
                records_regenerated.tsv; do
        cp "$PACKAGE/$name" "$HERE/" || true
    done
    cp "$PACKAGE"/tests/*.pytest.txt "$HERE/" || true
    cp "$PACKAGE"/tests/*.patch "$PACKAGE/records/INTEGRATION_RUN.md" "$HERE/" || true
    for name in README.md DEVIATIONS.md; do
        [ -s "$PACKAGE/$name" ] && cp "$PACKAGE/$name" "$HERE/$name.packaged"
    done
    # the tarball, unpacked elsewhere, against its own list of checksums
    local check; check=$(mktemp -d "$REGEN/unpack_check.XXXXXX")
    tar -xzf "$HERE/regenerated-$EXPECTED_VERSION.tar.gz" -C "$check"
    ( cd "$check/regenerated-$EXPECTED_VERSION"
      echo "PACKAGE_FILES $(find . -type f | wc -l)"
      echo "SHA256SUMS_ENTRIES $(wc -l < SHA256SUMS.txt)"
      echo "SHA256SUMS_FAILED $(sha256sum -c SHA256SUMS.txt 2>/dev/null | grep -vc ': OK$')"
      echo "UNPACKED_BYTES $(du -sb . | cut -f1)" )
    rm -rf "$check"
    echo "TARBALL_BYTES $(stat -c %s "$HERE/regenerated-$EXPECTED_VERSION.tar.gz")"
    echo "TARBALL_SHA256 $(sha256sum "$HERE/regenerated-$EXPECTED_VERSION.tar.gz" | cut -d' ' -f1)"
    [ -z "$(git -C "$REPO" status --porcelain)" ] && echo "WORKING_CLONE_CLEAN yes" || echo "WORKING_CLONE_CLEAN no"
}

chain_run() {   # two threads, then four, then six: never more than six at a time
    step input_checksums do_input_checksums
    step mock do_mock mock_results
    step mock_repeat do_mock mock_results_repeat
    step reference_results do_reference_results
    step low_copy do_low_copy
    step liba_gfa_4 do_liba_gfa 4
    step liba_gfa_6 do_liba_gfa 6
    step liba_verify do_liba_verify 4 6
    step consensus_nograph do_consensus_nograph
}

chain_package() {  # always made again from the files of the runs
    local name
    for name in reach evidence records tests integration compare finish summary manifest pack; do
        rm -f "$DONE/$name"
    done
    step reach do_reach
    step evidence do_evidence
    step records do_records
    step tests do_tests
    step integration do_integration
    step compare do_compare
    step finish do_finish
    step summary do_summary
    step manifest do_manifest
    step pack do_pack
}

if [ $# -eq 0 ]; then
    echo "usage: $0 prep | chain_run | chain_package | run <step function> [arguments]"
    exit 2
fi
for helper in "$HERE"/*.sh "$HERE"/*.py "$HERE"/*.md "$HERE"/*.tsv "$HERE"/*.patch; do
    [ -f "$helper" ] && [ "$HERE" != "$SCRIPTS" ] && cp "$helper" "$SCRIPTS/"
done
case "$1" in
    prep) step prep do_prep ;;
    chain_run) chain_run ;;
    chain_package) chain_package ;;
    run) shift; name=$1; shift; rm -f "$DONE/$name"; step "$name" "do_$name" "$@" ;;
    *) echo "unknown command $1"; exit 2 ;;
esac
