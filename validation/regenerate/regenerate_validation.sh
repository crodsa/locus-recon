#!/usr/bin/env bash
# Regenerate the validation analyses of Locus-Recon with the frozen release.
#
# Each analysis is one step.  A step writes its log to $REGEN/logs/<step>.log
# (with the commands, through `set -x`), its exit code and run time to
# $REGEN/logs/steps.tsv, and a marker to $REGEN/done/<step>; a step whose marker
# exists is skipped, so the script can be run again after an interruption.
# Thread counts and memory caps are those of the records of release 1.0.0,
# which the steps name, so that a difference between the two releases is not a
# difference of settings.
#
# Usage:
#   bash regenerate_validation.sh prep            # clone, check the digest, stage inputs
#   bash regenerate_validation.sh chain_small     # two threads
#   bash regenerate_validation.sh chain_main      # four threads, then six
#   bash regenerate_validation.sh chain_package   # block sizes, evidence, records, comparison, tests, package
#   bash regenerate_validation.sh run <step function> [arguments]   # one step
#
# Environment (defaults are the paths of the host used for the revision):
#   REGEN        working directory of the regeneration
#   SOURCE_REPO  git checkout of the release (cloned, never modified)
#   READ_CACHE   directory with ENA FASTQ files already downloaded (optional)
#   ANNOT_CACHE  directory with the RefSeq annotation FASTA files of the
#                1.0.0 run (optional; they are downloaded when absent)

REGEN=${REGEN:-/data/PROJECTS/locus_recon_revision/regen}
SOURCE_REPO=${SOURCE_REPO:-/data/PROJECTS/locus_recon_revision/repo}
READ_CACHE=${READ_CACHE:-/home/crodsa/claude_scratch/locus_recon_comparators/downloads}
ANNOT_CACHE=${ANNOT_CACHE:-/data/PROJECTS/locus-recon-validation/tier_work/downloads}
CONDA_SH=${CONDA_SH:-/home/crodsa/conda/etc/profile.d/conda.sh}
CONDA_ENV=${CONDA_ENV:-locus-recon-101}
COMPARATOR_ENV=${COMPARATOR_ENV:-lr-comparators}
EARLIER_TIER_WORK=${EARLIER_TIER_WORK:-/data/PROJECTS/locus-recon-validation/tier_work}
DATASETS_CLI=${DATASETS_CLI:-/home/crodsa/bin/datasets}
EXPECTED_DIGEST=${EXPECTED_DIGEST:-2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0}
EXPECTED_COMMIT=${EXPECTED_COMMIT:-884cdf3}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$REGEN/repo
WORK=$REGEN/work
LOGS=$REGEN/logs
DONE=$REGEN/done
INPUTS=$REGEN/inputs
RECORDS=$REGEN/records
SCRIPTS=$REGEN/scripts
PACKAGE=$REGEN/package/regenerated-1.0.1
PYDEPS=$REGEN/pydeps
HP_RUNS="SRR18335541 SRR18335539 SRR18335523 SRR18335501 SRR18335507"
LIBA_CONTIGS=".11940_5_90.1,.11940_5_90.29,.11940_5_90.56,.11940_5_90.59,.11940_5_90.63,.11940_5_90.76,.11940_5_90.77,.11940_5_90.8"

# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$CONDA_ENV"
set -uo pipefail
mkdir -p "$WORK" "$LOGS" "$DONE" "$INPUTS" "$RECORDS" "$SCRIPTS"
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

wait_for() {  # wait for the marker of another chain
    while [ ! -e "$DONE/$1" ]; do sleep 30; done
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

# ----------------------------------------------------------------------- prep
do_prep() {
    [ -d "$REPO/.git" ] || git clone -q "$SOURCE_REPO" "$REPO"
    cd "$REPO"
    git rev-parse HEAD
    git status --porcelain
    case "$(git rev-parse HEAD)" in
        "$EXPECTED_COMMIT"*) echo "COMMIT_OK" ;;
        *) echo "UNEXPECTED COMMIT"; exit 3 ;;
    esac
    local digest
    digest=$(python -c "from validation.common import source_sha256; print(source_sha256('.'))")
    echo "package digest of the working clone: $digest"
    [ "$digest" = "$EXPECTED_DIGEST" ] || { echo "DIGEST MISMATCH"; exit 3; }
    # the package installed in the environment must be the same code
    ( cd / && env -u PYTHONPATH python - <<'EOF'
import hashlib, pathlib, locus_recon
root = pathlib.Path(locus_recon.__file__).resolve().parent.parent
digest = hashlib.sha256()
for path in sorted((root / "locus_recon").rglob("*.py")):
    digest.update(path.relative_to(root).as_posix().encode() + b"\0")
    digest.update(path.read_bytes() + b"\0")
print("installed package:", locus_recon.__file__, locus_recon.VERSION, digest.hexdigest())
EOF
    )
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
        echo "pytest	$(python -c 'import pytest; print(pytest.__version__)' 2>&1 | tail -1)"
        echo "matplotlib	$(python -c 'import matplotlib; print(matplotlib.__version__)' 2>&1 | tail -1)"
        echo "git	$(git --version)"
        echo "platform	$(uname -srm)"
        echo "cpu	$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
    } > "$REGEN/TOOL_VERSIONS.tsv"
    cat "$REGEN/TOOL_VERSIONS.tsv"
    conda list --explicit > "$REGEN/conda_explicit_${CONDA_ENV}.txt" || true

    # inputs that are not downloaded by the runners
    [ -s "$INPUTS/LIBA6656_ST154.fasta" ] || cp "$HERE/LIBA6656_ST154.fasta" "$INPUTS/LIBA6656_ST154.fasta"
    cp "$INPUTS/LIBA6656_ST154.fasta" "$INPUTS/LIBA6656_ST154.assembly.fasta"
    link_reads "$INPUTS/reads" $HP_RUNS ERR467623
    ls -l "$INPUTS/reads"

    # fixed baits and genome of strain 26695 (NC_000915.1), through the function
    # of the low-copy runner; the completeness benchmark uses the same record
    python - <<EOF
from pathlib import Path
from validation.workflows.run_low_copy_23S import _prepare_fixed_baits
print(_prepare_fixed_baits(Path("$WORK/low_copy_23S")))
EOF

    # closed genome of MRSN 56 (GenBank GCA_019458485.1), reference of the aphA1 runs
    mkdir -p "$WORK/aphA1"
    if [ ! -s "$WORK/aphA1/MRSN56_GCA.fasta" ]; then
        ( cd "$WORK/aphA1" && "$DATASETS_CLI" download genome accession GCA_019458485.1 \
              --include genome --filename GCA_019458485.1.zip --no-progressbar \
          && unzip -o -q GCA_019458485.1.zip -d GCA_019458485.1_dataset \
          && cp GCA_019458485.1_dataset/ncbi_dataset/data/GCA_019458485.1/*_genomic.fna MRSN56_GCA.fasta )
    fi
    grep '>' "$WORK/aphA1/MRSN56_GCA.fasta"
    sha256sum "$WORK/aphA1/MRSN56_GCA.fasta" "$INPUTS/LIBA6656_ST154.fasta"

    # parameters of the 1.0.0 records that the steps mirror
    python - <<'EOF'
import json
from pathlib import Path
v = Path("validation")
for name in ("tier-calibration/PROVENANCE.json", "frame-and-graph-evidence/PROVENANCE.json"):
    record = json.loads((v / name).read_text())
    for key, block in (record.get("stages") or record.get("runs") or {}).items():
        print(name, key, json.dumps(block["parameters"], sort_keys=True))
for name in ("low-copy-23S/run_summary.json", "constructed-copy-series/run_summary.json",
             "aphA1-copy-number/run_summary.json", "liba6656-depth/run_summary.json"):
    print(name, json.dumps(json.loads((v / name).read_text())["parameters"], sort_keys=True))
record = json.loads((v / "liba6656-gfa-rerun/RELEASE_VERIFICATION.json").read_text())
print("liba6656-gfa-rerun", [run["parameters"] for run in record["runs"]])
EOF
    sed -n 1,80p validation/depth-copy-number/overlap_null/README.md
}

do_input_checksums() {
    cd "$INPUTS"
    : > INPUT_SHA256.tsv
    local path
    for path in reads/*.fastq.gz LIBA6656_ST154.fasta; do
        printf '%s\t%s\t%s\n' "$path" "$(stat -L -c %s "$path")" "$(sha256sum "$path" | cut -d' ' -f1)" >> INPUT_SHA256.tsv
    done
    cat INPUT_SHA256.tsv
}

# ----------------------------------------------------------------------- mock
do_mock_dataset_check() {
    cd "$REPO"
    rm -rf "$WORK/mock/regenerated_dataset"
    mkdir -p "$WORK/mock"
    python - <<EOF
import sys
sys.path.insert(0, "validation")
from pathlib import Path
from generate_mock_dataset import build_dataset
build_dataset(Path("$WORK/mock/regenerated_dataset"), force=False)
EOF
    if diff -r validation/mock_dataset "$WORK/mock/regenerated_dataset"; then
        echo "MOCK_DATASET_IDENTICAL_TO_COMMITTED yes"
    else
        echo "MOCK_DATASET_IDENTICAL_TO_COMMITTED no"
    fi
}

do_mock() {  # do_mock <results directory name>
    cd "$REPO"
    python validation/run_mock_validation.py \
        --results-dir "$WORK/mock/$1" --threads 2 --memory-gb 4 --force
}

# --------------------------------------------------------------- completeness
do_completeness() {
    cd "$REPO"
    local baits=$WORK/low_copy_23S/fixed_baits
    rm -rf "$WORK/completeness"
    python validation/run_completeness_benchmark.py \
        --bait "$baits/NC_000915.1_23S_gyrB.fasta" --locus-id gyrB \
        --genome "$baits/NC_000915.1.gbff" \
        --outdir "$WORK/completeness/from_genome"
    python validation/run_completeness_benchmark.py \
        --from-deposit validation/completeness-benchmark \
        --outdir "$WORK/completeness/from_deposit"
    # are the constructed assemblies those of the deposit?
    local case same=0 other=0
    for case in validation/completeness-benchmark/cases/*; do
        if cmp -s "$case/assembly.fasta" "$WORK/completeness/from_genome/cases/$(basename "$case")/assembly.fasta"; then
            same=$((same + 1))
        else
            other=$((other + 1)); echo "assembly differs: $(basename "$case")"
        fi
    done
    echo "CONSTRUCTED_ASSEMBLIES identical=$same different=$other"
}

# The records of the two benchmarks replace those of 1.0.0 in the working
# clone, because stage A of the tier calibration reads them from there.
do_overlay_benchmarks() {
    cd "$REPO"
    python "$SCRIPTS/make_reference_results.py" \
        --results-dir "$WORK/mock/mock_results" --out "$RECORDS/reference_results"
    cp "$RECORDS/reference_results/locus_recon_report_mockLocus.tsv" \
       "$RECORDS/reference_results/validation_results.tsv" \
       "$RECORDS/reference_results/validation_summary.json" validation/reference_results/
    rm -rf "$RECORDS/completeness-benchmark"
    mkdir -p "$RECORDS/completeness-benchmark"
    cp -r "$WORK/completeness/from_genome/." "$RECORDS/completeness-benchmark/"
    cp "$RECORDS/completeness-benchmark/completeness_benchmark_results.tsv" \
       "$RECORDS/completeness-benchmark/completeness_benchmark_results.json" \
       validation/completeness-benchmark/
    git status --porcelain
}

# --------------------------------------------------------------- low-copy 23S
do_low_copy() {
    cd "$REPO"
    link_reads "$WORK/low_copy_23S/reads" $HP_RUNS
    # As in the 1.0.0 record (its run_summary.json, execution_note), the five
    # SPAdes --isolate drafts are those of the closed-genome stage of the tier
    # calibration: the runner keeps a draft that already exists.
    local strain
    for strain in Hpfe0001 Hpfe0002 Hpfe0003 Hpfe0004 Hpfe0006; do
        rm -rf "$WORK/low_copy_23S/samples/$strain"
        mkdir -p "$WORK/low_copy_23S/samples/$strain/spades"
        cp "$WORK/tier/$strain/spades/scaffolds.fasta" "$WORK/tier/$strain/spades/contigs.fasta" \
           "$WORK/tier/$strain/spades/assembly_graph_after_simplification.gfa" \
           "$WORK/tier/$strain/spades/params.txt" "$WORK/tier/$strain/spades/spades.log" \
           "$WORK/low_copy_23S/samples/$strain/spades/"
    done
    set +e
    python validation/workflows/run_low_copy_23S.py \
        --work-dir "$WORK/low_copy_23S" --threads 4 --memory-gb 8
    local rc=$?
    set -e
    # exit 1 is the outcome when a primary criterion is not met; evaluation.json
    # names the criterion (as in 1.0.0: the two-copy consensus of Hpfe0004 is a
    # lower bound)
    echo "RUN_LOW_COPY_23S_EXIT $rc"
}

# ----------------------------------------------------------- tier calibration
do_tier() {
    cd "$REPO"
    mkdir -p "$WORK/tier/downloads"
    link_reads "$WORK/tier/downloads" $HP_RUNS
    if [ -d "$ANNOT_CACHE" ]; then
        cp -n "$ANNOT_CACHE"/*_from_genomic.fna.gz "$WORK/tier/downloads/" || true
    fi
    ls -l "$WORK/tier/downloads"
    # The runner keeps the tables of stages it does not run.  All stages are
    # run here, from a folder without the result files of 1.0.0.
    ( cd validation/tier-calibration && rm -f TIER_CALIBRATION.json PROVENANCE.json \
        INPUT_CHECKSUMS.tsv tier_calibration*.tsv tier_calibration_tables.md \
        tier_calibration.pdf tier_calibration.svg bait_*_26695.fasta locus_selection.tsv )
    python validation/run_tier_calibration.py \
        --outdir validation/tier-calibration \
        --manifest validation/tier-calibration/closed_genome_manifest.tsv \
        --workdir "$WORK/tier" --threads 4 --memory-gb 14 \
        --divergent-loci 3 --depth-series 15,8 --depth-isolate
    rm -rf "$RECORDS/tier-calibration"
    cp -r validation/tier-calibration "$RECORDS/tier-calibration"
}

# ---------------------------------------------------------------------- aphA1
do_aphA1_download() {
    cd "$REPO"
    python - <<EOF
from pathlib import Path
from validation.workflows.run_aphA1_validation import _download_ena_pair, _load_manifest
for row in _load_manifest(Path("validation/manifests/aphA1_copy_number.tsv")):
    for record in _download_ena_pair(row["illumina_run"], Path("$WORK/aphA1/reads")):
        print(record, flush=True)
EOF
}

# The downloader of the runner keeps a file that a dropped connection left
# incomplete, and then stops at its MD5.  This step fetches the same ENA files
# with resumption and checks each against the MD5 that ENA publishes; the
# runner verifies them again and writes download_provenance.json.
fetch_ena_file() {  # fetch_ena_file <name> <url> <md5> <bytes>
    local name=$1 url=$2 md5=$3 bytes=$4
    local dest=$WORK/aphA1/reads/$name attempt size
    if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then
        echo "present and verified: $name"
        return 0
    fi
    rm -f "$dest"
    for attempt in 1 2 3 4 5 6 7 8 9 10 11 12; do
        curl -sS -L --fail -C - --connect-timeout 60 --speed-time 120 --speed-limit 1000 \
            -o "$dest.part" "$url" || true
        size=$(stat -c %s "$dest.part" 2>/dev/null || echo 0)
        if [ "$size" = "$bytes" ]; then
            if [ "$(md5sum "$dest.part" | cut -d' ' -f1)" = "$md5" ]; then
                mv "$dest.part" "$dest"
                echo "downloaded and verified: $name ($bytes bytes, attempt $attempt)"
                return 0
            fi
            echo "MD5 mismatch at full size for $name; starting again"
            rm -f "$dest.part"
        elif [ "$size" -gt "$bytes" ]; then
            rm -f "$dest.part"
        else
            echo "attempt $attempt: $name has $size of $bytes bytes"
        fi
        sleep 15
    done
    echo "FAILED: $name"
    return 1
}

do_aphA1_fetch() {  # do_aphA1_fetch <downloads at a time>
    cd "$REPO"
    mkdir -p "$WORK/aphA1/reads"
    rm -f "$WORK/aphA1/reads"/*.fastq.fastq.gz.part
    python - > "$WORK/aphA1/ena_files.tsv" <<'EOF'
import urllib.parse
import urllib.request
from pathlib import Path
from validation.workflows.run_aphA1_validation import ENA_REPORT, _load_manifest
for row in _load_manifest(Path("validation/manifests/aphA1_copy_number.tsv")):
    query = urllib.parse.urlencode({
        "accession": row["illumina_run"], "result": "read_run",
        "fields": "fastq_ftp,fastq_md5,fastq_bytes", "format": "tsv"})
    with urllib.request.urlopen(f"{ENA_REPORT}?{query}", timeout=120) as response:
        lines = response.read().decode().strip().splitlines()
    report = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    files = zip(report["fastq_ftp"].split(";"), report["fastq_md5"].split(";"),
                report["fastq_bytes"].split(";"))
    for mate, (url, md5, size) in enumerate(files, start=1):
        print(f"{row['illumina_run']}_{mate}.fastq.gz\thttps://{url}\t{md5}\t{size}")
EOF
    cat "$WORK/aphA1/ena_files.tsv"
    local name url md5 bytes running=0 failed=0
    while IFS=$'\t' read -r name url md5 bytes; do
        fetch_ena_file "$name" "$url" "$md5" "$bytes" &
        running=$((running + 1))
        if [ "$running" -ge "${1:-3}" ]; then
            wait -n || failed=$((failed + 1))
            running=$((running - 1))
        fi
    done < "$WORK/aphA1/ena_files.tsv"
    while [ "$running" -gt 0 ]; do
        wait -n || failed=$((failed + 1))
        running=$((running - 1))
    done
    ls -l "$WORK/aphA1/reads"
    echo "APHA1_FETCH failed=$failed"
    [ "$failed" -eq 0 ]
}

do_aphA1() {
    cd "$REPO"
    set +e
    python validation/workflows/run_aphA1_validation.py \
        --reference "$WORK/aphA1/MRSN56_GCA.fasta" \
        --reads-dir "$WORK/aphA1/reads" \
        --output-dir "$WORK/aphA1/results" --threads 4 --download-missing
    local rc=$?
    set -e
    echo "RUN_APHA1_EXIT $rc"
}

# ------------------------------------------------------ constructed copy series
do_constructed() {
    cd "$REPO"
    python validation/workflows/run_constructed_copy_series.py \
        --work-dir "$WORK/constructed_copy_series" --threads 4 \
        --deposit-dir "$RECORDS/constructed-copy-series"
    cp validation/constructed-copy-series/PRESPECIFICATION.md "$RECORDS/constructed-copy-series/"
}

# ------------------------------------------------------------------ LIBA-6656
do_liba_depth() {
    cd "$REPO"
    link_reads "$WORK/liba6656_depth/reads" ERR467623
    python validation/workflows/run_liba6656_depth.py \
        --draft-assembly "$INPUTS/LIBA6656_ST154.assembly.fasta" \
        --reads-dir "$WORK/liba6656_depth/reads" \
        --output-dir "$WORK/liba6656_depth/results" --threads 4
}

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
}

do_overlap_null() {
    cd "$REPO"
    python validation/depth-copy-number/overlap_null/run_overlap_null.py \
        --assembly "$INPUTS/LIBA6656_ST154.assembly.fasta" \
        --bam "$WORK/liba6656_depth/results/ERR467623_vs_LIBA6656_draft.sorted.bam" \
        --locus-contigs "$LIBA_CONTIGS" \
        --case-evaluation "$WORK/liba6656_depth/results/evaluation.json" \
        --outdir "$WORK/overlap_null" --threads 4
}

# ------------------------------------------------------ frame and graph checks
do_frame_graph() {  # do_frame_graph <threads>
    cd "$REPO"
    rm -rf "$WORK/frame_graph"
    mkdir -p "$WORK/frame_graph"
    python validation/run_frame_and_graph_validation.py \
        --deposit validation/liba6656-gfa-rerun \
        --outdir "$WORK/frame_graph/frame-and-graph-evidence" \
        --reads-r1 "$INPUTS/reads/ERR467623_1.fastq.gz" \
        --reads-r2 "$INPUTS/reads/ERR467623_2.fastq.gz" \
        --threads "$1"
}

# ----------------------------------------------- block size of the bootstrap
do_block_size() {  # do_block_size <cases evaluated at a time>
    cd "$REPO"
    [ -s validation/run_block_size_sensitivity.py ] \
        || cp "$SCRIPTS/run_block_size_sensitivity.py" validation/run_block_size_sensitivity.py
    python validation/run_block_size_sensitivity.py \
        --aphA1-results "$WORK/aphA1/results" \
        --aphA1-reference "$WORK/aphA1/MRSN56_GCA.fasta" \
        --low-copy-work "$WORK/low_copy_23S" \
        --output "$WORK/block_size/block_size_sensitivity.tsv" --jobs "${1:-4}"
}

# ------------------------------------- packages that the environment lacks
# pytest and matplotlib are installed beside the environment, not into it.
do_pydeps() {
    python -m pip install --quiet --target "$PYDEPS/pytest" pytest
    python -m pip install --quiet --target "$PYDEPS/matplotlib" matplotlib
    PYTHONPATH=$PYDEPS/pytest python -c "import pytest; print('pytest', pytest.__version__)"
    PYTHONPATH=$PYDEPS/matplotlib python -c "import matplotlib; print('matplotlib', matplotlib.__version__)"
}

do_tier_figure() {
    cd "$REPO"
    PYTHONPATH=$PYDEPS/matplotlib:$REPO python - <<EOF
import sys
sys.path.insert(0, "validation")
import matplotlib
import run_tier_calibration as calibration
combined, missing, rows = calibration.combine("$RECORDS/tier-calibration")
print("missing stages:", missing)
print("figure:", calibration.plot_calibration(rows, "$RECORDS/tier-calibration"))
print("matplotlib", matplotlib.__version__)
EOF
}

# -------------------------------------------------------------------- package
do_evidence() {
    python "$SCRIPTS/collect_regenerated.py" evidence --regen "$REGEN" --package "$PACKAGE" \
        --repo "$REPO" --earlier-tier-work "$EARLIER_TIER_WORK"
}

do_comparator() {
    cd "$REPO"
    if ! grep -q -- "--locus-recon-evidence" validation/run_comparator_benchmark.py; then
        python "$SCRIPTS/patch_comparator.py" \
            --source "$SOURCE_REPO/validation/run_comparator_benchmark.py" \
            --output validation/run_comparator_benchmark.py
    fi
    cp validation/run_comparator_benchmark.py "$SCRIPTS/run_comparator_benchmark.py"
    local work=$WORK/comparator name
    rm -rf "$work"
    mkdir -p "$work/downloads"
    # ARIBA and KMA are not run again: their deposited raw outputs are scored
    tar -xzf "$SOURCE_REPO/validation/comparator-benchmark/raw_tool_outputs.tar.gz" -C "$work"
    find "$work" -maxdepth 3 | sort
    cp "$WORK/tier/downloads"/*_from_genomic.fna.gz "$work/downloads/"
    set +u
    conda activate "$COMPARATOR_ENV"
    set -u
    python --version
    blastn -version
    # control: the released script, on the tables of 1.0.0, must reproduce the deposit
    python "$SOURCE_REPO/validation/run_comparator_benchmark.py" \
        --workdir "$work" --outdir "$work/control_released_script" --score-only
    for name in comparator_benchmark_results.tsv comparator_benchmark_summary.tsv \
                comparator_returned_sequences.fasta comparator_truth_sequences.fasta; do
        if cmp -s "$work/control_released_script/$name" \
                  "$SOURCE_REPO/validation/comparator-benchmark/$name"; then
            echo "CONTROL $name identical to the deposit of 1.0.0"
        else
            echo "CONTROL $name DIFFERS from the deposit of 1.0.0"
        fi
    done
    rm -rf "$PACKAGE/records/comparator-benchmark"
    python validation/run_comparator_benchmark.py \
        --workdir "$work" --outdir "$PACKAGE/records/comparator-benchmark" --score-only \
        --locus-recon-evidence "$PACKAGE/evidence"
}

do_records() {
    python "$SCRIPTS/collect_regenerated.py" records --regen "$REGEN" --package "$PACKAGE" \
        --old "$SOURCE_REPO/validation" --repo "$REPO" \
        --liba-readme-section "$SCRIPTS/liba6656_readme_section.md"
    mkdir -p "$PACKAGE/records/block-size-sensitivity"
    cp "$WORK/block_size/block_size_sensitivity.tsv" "$PACKAGE/records/block-size-sensitivity/"
}

do_compare() {
    python "$SCRIPTS/collect_regenerated.py" compare --regen "$REGEN" --package "$PACKAGE" \
        --old "$SOURCE_REPO/validation" --repo "$REPO" --notes "$SCRIPTS/difference_notes.tsv"
}

# The test suite of the repository, on three trees: the release as committed,
# the release with the regenerated LIBA-6656 record, and the release with
# every regenerated record copied over validation/.
do_tests() {
    local base=$WORK/tests tree
    rm -rf "$base"
    mkdir -p "$base"
    for tree in release_as_committed release_with_liba6656_record release_with_all_records \
                release_with_all_records_and_test_update; do
        git clone -q "$SOURCE_REPO" "$base/$tree"
        case "$tree" in
            release_with_liba6656_record)
                cp -r "$PACKAGE/records/liba6656-gfa-rerun/." "$base/$tree/validation/liba6656-gfa-rerun/" ;;
            release_with_all_records)
                cp -r "$PACKAGE/records/." "$base/$tree/validation/" ;;
            release_with_all_records_and_test_update)
                cp -r "$PACKAGE/records/." "$base/$tree/validation/"
                # one test pins the SHA-256 of the three aphA1 record files of
                # 1.0.0; the suggested update replaces a pinned value only
                # where the regenerated file differs
                ( cd "$base/$tree" && python - <<'EOF'
import hashlib, subprocess
from pathlib import Path
test = Path("tests/test_validation_workflows.py")
text = test.read_text()
for name in ("results.tsv", "evaluation.json", "run_summary.json"):
    path = f"validation/aphA1-copy-number/{name}"
    old = hashlib.sha256(subprocess.run(["git", "show", f"HEAD:{path}"], capture_output=True,
                                        check=True).stdout).hexdigest()
    new = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    print(name, old, "->", new, "(unchanged)" if old == new else "")
    if old != new:
        assert text.count(old) == 1, name
        text = text.replace(old, new)
test.write_text(text)
EOF
                  git diff -- tests > "$base/suggested_test_update.patch" ) ;;
        esac
        ( cd "$base/$tree"
          git status --porcelain | wc -l
          set +e
          env PYTHONPATH="$PYDEPS/pytest" python -m pytest -q -p no:cacheprovider \
              > "$base/$tree.pytest.txt" 2>&1
          echo "PYTEST_EXIT $tree $?"
          set -e
          tail -n 25 "$base/$tree.pytest.txt" )
    done
}

do_manifest() {
    mkdir -p "$PACKAGE/scripts/validation/regenerate" "$PACKAGE/logs" "$PACKAGE/inputs/refseq_annotation" "$PACKAGE/tests"
    cp "$SCRIPTS/regenerate_validation.sh" "$SCRIPTS/collect_regenerated.py" \
       "$SCRIPTS/make_reference_results.py" "$SCRIPTS/patch_comparator.py" \
       "$PACKAGE/scripts/validation/regenerate/"
    cp "$SCRIPTS/run_block_size_sensitivity.py" "$SCRIPTS/run_comparator_benchmark.py" \
       "$PACKAGE/scripts/validation/"
    find "$LOGS" -maxdepth 1 -type f -size -5000k -exec cp {} "$PACKAGE/logs/" \;
    cp "$INPUTS/INPUT_SHA256.tsv" "$REGEN/TOOL_VERSIONS.tsv" "$REGEN"/conda_explicit_*.txt "$PACKAGE/inputs/"
    cp "$WORK/tier/downloads"/*_from_genomic.fna.gz "$PACKAGE/inputs/refseq_annotation/"
    cp "$WORK/tests"/*.pytest.txt "$WORK/tests/suggested_test_update.patch" "$PACKAGE/tests/" || true
    if [ -s "$SCRIPTS/readme_edits.py" ]; then
        cp "$SCRIPTS/readme_edits.py" "$PACKAGE/scripts/validation/regenerate/"
    fi
    if [ -s "$SCRIPTS/difference_notes.tsv" ]; then
        cp "$SCRIPTS/difference_notes.tsv" "$SCRIPTS/liba6656_readme_section.md" \
           "$PACKAGE/scripts/validation/regenerate/"
    fi
    local name
    for name in README.md DEVIATIONS.md; do
        if [ -s "$SCRIPTS/$name" ]; then cp "$SCRIPTS/$name" "$PACKAGE/$name"; fi
    done
    cp "$PACKAGE/records/block-size-sensitivity/block_size_sensitivity.tsv" "$PACKAGE/block_size_sensitivity.tsv"
    python "$SCRIPTS/collect_regenerated.py" summary --regen "$REGEN" --package "$PACKAGE" \
        --repo "$REPO" --digest "$EXPECTED_DIGEST" --commit "$(git -C "$REPO" rev-parse HEAD)"
    python "$SCRIPTS/collect_regenerated.py" integration --regen "$REGEN" --package "$PACKAGE" \
        --repo "$REPO" --digest "$EXPECTED_DIGEST" --commit "$(git -C "$REPO" rev-parse HEAD)"
    cp "$PACKAGE/records/block-size-sensitivity/block_size_sensitivity.tsv" "$PACKAGE/block_size_sensitivity.tsv"
    python "$SCRIPTS/collect_regenerated.py" manifest --regen "$REGEN" --package "$PACKAGE" \
        --repo "$REPO" --digest "$EXPECTED_DIGEST" --commit "$(git -C "$REPO" rev-parse HEAD)"
}

do_pack() {
    cd "$REGEN/package"
    tar -czf "$HERE/regenerated-validation.tar.gz" regenerated-1.0.1
    cp "$PACKAGE/old_vs_new.tsv" "$PACKAGE/old_vs_new_summary.tsv" "$PACKAGE/block_size_sensitivity.tsv" \
       "$PACKAGE/PROVENANCE.json" "$PACKAGE/SHA256SUMS.txt" "$PACKAGE/FILES_LEFT_ON_HOST.tsv" \
       "$PACKAGE/exactness_check.tsv" "$PACKAGE/tests/suggested_test_update.patch" "$HERE/" || true
    cp "$SCRIPTS/run_comparator_benchmark.py" "$PACKAGE/records/INTEGRATION_RUN.md" \
       "$PACKAGE/records/reference_results/BENCHMARK_REPRODUCTION.md" "$HERE/" || true
    cp "$PACKAGE/summary.tsv" "$HERE/" || true
    ls -l "$HERE"
    tar -tzf "$HERE/regenerated-validation.tar.gz" | wc -l
}

# Check 6 maps the reads to 512 nearly identical paths.  Its record of 1.0.0
# was made with 24 threads.  The check is repeated with 2, 4 and 6 threads, and
# every table of path scores is compared with every other and with that record.
do_frame_graph_threads() {
    cd "$REPO"
    local threads
    for threads in 2 4 6; do
        rm -rf "$WORK/frame_graph_threads/threads_$threads"
        mkdir -p "$WORK/frame_graph_threads"
        python validation/run_frame_and_graph_validation.py \
            --deposit validation/liba6656-gfa-rerun \
            --outdir "$WORK/frame_graph_threads/threads_$threads" \
            --reads-r1 "$INPUTS/reads/ERR467623_1.fastq.gz" \
            --reads-r2 "$INPUTS/reads/ERR467623_2.fastq.gz" \
            --threads "$threads" > "$WORK/frame_graph_threads/threads_$threads.log" 2>&1
    done
    python - <<EOF
import csv, gzip, hashlib, itertools
from pathlib import Path
work = Path("$WORK")
runs = {
    "release 1.0.0 record (24 threads)": Path("$SOURCE_REPO/validation/frame-and-graph-evidence"),
    "regenerated record (6 threads)": work / "frame_graph" / "frame-and-graph-evidence",
    "repetition with 2 threads": work / "frame_graph_threads" / "threads_2",
    "repetition with 4 threads": work / "frame_graph_threads" / "threads_4",
    "repetition with 6 threads": work / "frame_graph_threads" / "threads_6",
}
def table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
rows = []
for (first, a), (second, b) in itertools.combinations(runs.items(), 2):
    old, new = table(a / "graph_path_scores.tsv"), table(b / "graph_path_scores.tsv")
    differing, columns, largest = 0, set(), 0
    for x, y in zip(old, new):
        changed = [key for key in x if x[key] != y.get(key)]
        if changed:
            differing += 1
            columns.update(changed)
            if "mapped_reads" in changed:
                largest = max(largest, abs(int(x["mapped_reads"]) - int(y["mapped_reads"])))
    rows.append(["graph_path_scores.tsv", first, second, len(old), differing,
                 ";".join(sorted(columns)), largest])
    same = (a / "VALIDATION_SUMMARY.json").read_bytes() == (b / "VALIDATION_SUMMARY.json").read_bytes()
    rows.append(["VALIDATION_SUMMARY.json", first, second, 1, 0 if same else 1, "", ""])
    for mate in ("R1", "R2"):
        name = f"scoring_positive_control_{mate}.fastq.gz"
        digests = [hashlib.sha256(gzip.open(folder / name).read()).hexdigest() for folder in (a, b)]
        rows.append([f"{name} (reads inside the gzip file)", first, second, 1,
                     0 if digests[0] == digests[1] else 1, "", ""])
with open(work / "frame_graph" / "thread_count_comparison.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["file", "run_a", "run_b", "rows_compared", "rows_differing",
                     "columns_differing", "largest_difference_in_mapped_reads"])
    writer.writerows(rows)
for row in rows:
    if row[0] == "graph_path_scores.tsv":
        print(row)
EOF
}

# --------------------------------------------------------------------- chains
chain_small() {   # two threads at most
    step input_checksums do_input_checksums
    step mock_dataset_check do_mock_dataset_check
    step mock do_mock mock_results
    step mock_repeat do_mock mock_results_repeat
    step completeness do_completeness
    step overlay_benchmarks do_overlay_benchmarks
    step aphA1_download do_aphA1_download
    echo "chain_small finished $(date -u +%FT%TZ)" > "$DONE/chain_small_finished"
}

chain_main() {    # four threads while chain_small may run, six afterwards
    wait_for overlay_benchmarks
    step tier do_tier
    step low_copy do_low_copy
    wait_for aphA1_download
    step aphA1 do_aphA1
    step constructed do_constructed
    step liba_depth do_liba_depth
    step liba_gfa_4 do_liba_gfa 4
    step overlap_null do_overlap_null
    wait_for chain_small_finished
    step liba_gfa_6 do_liba_gfa 6
    step liba_verify do_liba_verify 4 6
    step frame_graph do_frame_graph 6
    echo "chain_main finished $(date -u +%FT%TZ)" > "$DONE/chain_main_finished"
}

chain_package() { # after the two chains; four processes at most
    # the assembly of the package is always done again; the steps before it are kept
    for name in evidence comparator records compare tests manifest pack; do
        rm -f "$DONE/$name"
    done
    step pydeps do_pydeps
    step tier_figure do_tier_figure
    step block_size do_block_size 4
    step frame_graph_threads do_frame_graph_threads
    step evidence do_evidence
    step comparator do_comparator
    step records do_records
    step compare do_compare
    step tests do_tests
    step manifest do_manifest
    step pack do_pack
}

if [ $# -eq 0 ]; then
    echo "usage: $0 prep | chain_small | chain_main | chain_package | run <step function> [arguments]"
    exit 2
fi
for helper in "$HERE"/*.sh "$HERE"/*.py "$HERE"/*.md "$HERE"/*.tsv; do
    [ -f "$helper" ] && [ "$HERE" != "$SCRIPTS" ] && cp "$helper" "$SCRIPTS/"
done
case "$1" in
    prep) step prep do_prep ;;
    chain_small|chain_main|chain_package) "$1" ;;
    run) shift; name=$1; shift; step "${STEP_NAME:-${name#do_}}" "$name" "$@" ;;
    redo) shift; name=$1; shift; rm -f "$DONE/${STEP_NAME:-${name#do_}}"
          step "${STEP_NAME:-${name#do_}}" "$name" "$@" ;;
    *) echo "unknown command: $1"; exit 2 ;;
esac
ls -l "$DONE"
tail -n 40 "$LOGS/steps.tsv" 2>/dev/null || true
