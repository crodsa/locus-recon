#!/usr/bin/env python3
"""Run the prespecified external two-copy 23S rRNA validation."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Bio import SeqIO

from locus_recon.assembly_graph import parse_gfa
from locus_recon.copy_number import (
    discover_graph_target_hits,
    infer_graph_contexts,
    reconcile_copy_number,
)
from locus_recon.depth_ratio import estimate_locus_copy_number
from validation.common import build_provenance, require_files, sha256_file, write_json


LEGACY_RESULT_FIELDS = (
    "strain", "biosample", "illumina_run", "assembly_accession", "locus",
    "truth_copies", "assembly_region_count", "ratio", "ratio_ci_low",
    "ratio_ci_high", "depth_call", "dosage_estimate", "dosage_status",
    "dosage_flags", "absolute_error", "ambiguity_index",
    "gene_total_aligned_bp", "gene_union_covered_bp",
    "gene_overlap_bp", "locus_regions", "region_gene_spans",
)


FIXED_REFERENCE = "NC_000915.1"
FIXED_REFERENCE_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=NC_000915.1&rettype=gbwithparts&retmode=text"
)
ENA_REPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "locus-recon-validation/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    if temporary.stat().st_size == 0:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"download was empty: {url}")
    temporary.replace(destination)


def _md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _ena_fastqs(run: str) -> list[tuple[str, str]]:
    query = urllib.parse.urlencode({
        "accession": run,
        "result": "read_run",
        "fields": "fastq_ftp,fastq_md5",
        "format": "tsv",
    })
    with urllib.request.urlopen(f"{ENA_REPORT}?{query}", timeout=60) as response:
        lines = response.read().decode().strip().splitlines()
    if len(lines) != 2:
        raise RuntimeError(f"ENA returned {len(lines) - 1} records for {run}, expected one")
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    urls, md5s = row["fastq_ftp"].split(";"), row["fastq_md5"].split(";")
    if len(urls) != 2 or len(md5s) != 2:
        raise RuntimeError(f"{run} is not represented by exactly two ENA FASTQ files")
    return [(f"https://{url}", md5) for url, md5 in zip(urls, md5s)]


def _prepare_reads(run: str, reads_dir: Path) -> tuple[Path, Path, list[dict]]:
    provenance: list[dict] = []
    paths: list[Path] = []
    for mate, (url, expected_md5) in enumerate(_ena_fastqs(run), start=1):
        destination = reads_dir / f"{run}_{mate}.fastq.gz"
        if not destination.exists():
            _download(url, destination)
        observed_md5 = _md5(destination)
        if observed_md5 != expected_md5:
            raise RuntimeError(
                f"MD5 mismatch for {destination}: expected {expected_md5}, observed {observed_md5}"
            )
        paths.append(destination)
        provenance.append({
            "run": run,
            "mate": mate,
            "url": url,
            "ena_md5": expected_md5,
            "local_md5": observed_md5,
            "sha256": sha256_file(destination),
            "size_bytes": destination.stat().st_size,
        })
    return paths[0], paths[1], provenance


def _prepare_fixed_baits(work_dir: Path) -> tuple[Path, dict[str, int]]:
    gbff = work_dir / "fixed_baits" / f"{FIXED_REFERENCE}.gbff"
    if not gbff.exists():
        _download(FIXED_REFERENCE_URL, gbff)
    record = SeqIO.read(gbff, "genbank")
    selected = {}
    for feature in record.features:
        product = " ".join(feature.qualifiers.get("product", []))
        gene = feature.qualifiers.get("gene", [""])[0]
        if feature.type == "rRNA" and "23S ribosomal RNA" in product and "23S" not in selected:
            selected["23S"] = feature.extract(record.seq)
        if feature.type == "CDS" and gene == "gyrB" and "gyrB" not in selected:
            selected["gyrB"] = feature.extract(record.seq)
    if set(selected) != {"23S", "gyrB"}:
        raise RuntimeError(f"could not extract both fixed baits from {FIXED_REFERENCE}")
    bait_path = work_dir / "fixed_baits" / f"{FIXED_REFERENCE}_23S_gyrB.fasta"
    with bait_path.open("w") as handle:
        for name in ("23S", "gyrB"):
            handle.write(f">{name}|source={FIXED_REFERENCE}\n")
            sequence = str(selected[name])
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    return bait_path, {name: len(sequence) for name, sequence in selected.items()}


def _audit_truth(row: dict[str, str], truth_dir: Path) -> tuple[dict, Path]:
    prefix = row["assembly_ftp"].rstrip("/")
    stem = prefix.rsplit("/", 1)[1]
    gbff = truth_dir / f"{stem}_genomic.gbff.gz"
    if not gbff.exists():
        _download(f"{prefix}/{stem}_genomic.gbff.gz", gbff)
    with gzip.open(gbff, "rt") as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    def is_plasmid(record) -> bool:
        source_features = [feature for feature in record.features if feature.type == "source"]
        has_plasmid_qualifier = any(
            feature.qualifiers.get("plasmid") for feature in source_features
        )
        description = str(record.description).lower()
        return has_plasmid_qualifier or " plasmid " in f" {description} "

    chromosomes = [
        record for record in records
        if record.annotations.get("topology") == "circular" and not is_plasmid(record)
    ]
    if len(chromosomes) != 1:
        raise RuntimeError(
            f"{row['assembly_accession']}: expected one circular non-plasmid chromosome, "
            f"observed {len(chromosomes)} among {len(records)} deposited replicons"
        )
    chromosome = chromosomes[0]
    other_replicons = [
        {
            "accession": record.id,
            "length": len(record),
            "topology": record.annotations.get("topology"),
            "classified_as_plasmid": is_plasmid(record),
        }
        for record in records if record is not chromosome
    ]
    rrna = []
    gyrb = []
    for feature in chromosome.features:
        product = " ".join(feature.qualifiers.get("product", []))
        gene = feature.qualifiers.get("gene", [""])[0]
        if feature.type == "rRNA" and "23S ribosomal RNA" in product:
            rrna.append(feature)
        if feature.type == "CDS" and gene == "gyrB":
            gyrb.append(feature)
    full_rrna = [feature for feature in rrna if len(feature) >= 2800]
    full_gyrb = [feature for feature in gyrb if len(feature) >= 2200]
    if len(full_rrna) != int(row["truth_23S_copies"]):
        raise RuntimeError(
            f"{row['assembly_accession']}: observed {len(full_rrna)} full 23S features"
        )
    if len(full_gyrb) != int(row["truth_control_copies"]):
        raise RuntimeError(
            f"{row['assembly_accession']}: observed {len(full_gyrb)} full gyrB CDS features"
        )
    return ({
        "strain": row["strain"],
        "biosample": row["biosample"],
        "assembly_accession": row["assembly_accession"],
        "chromosome_accession": chromosome.id,
        "chromosome_length": len(chromosome),
        "chromosome_topology": chromosome.annotations.get("topology"),
        "full_length_23S_copies": len(full_rrna),
        "23S_lengths": [len(feature) for feature in full_rrna],
        "full_length_gyrB_copies": len(full_gyrb),
        "gyrB_lengths": [len(feature) for feature in full_gyrb],
        "other_deposited_replicons": other_replicons,
        "gbff_sha256": sha256_file(gbff),
    }, gbff)


def _assemble_draft(
    r1: Path, r2: Path, output_dir: Path, *, threads: int, memory_gb: int, spades: str,
) -> Path:
    assembly = output_dir / "scaffolds.fasta"
    if assembly.exists() and assembly.stat().st_size > 0:
        return assembly
    if (output_dir / "params.txt").exists():
        command = [spades, "--continue", "-o", str(output_dir)]
    else:
        command = [
            spades, "--isolate", "-1", str(r1), "-2", str(r2), "-o", str(output_dir),
            "--threads", str(threads), "--memory", str(memory_gb),
        ]
    completed = subprocess.run(command, capture_output=True, text=True)
    (output_dir.parent / "spades.command.txt").write_text(" ".join(command) + "\n")
    (output_dir.parent / "spades.stderr.log").write_text(completed.stderr)
    if completed.returncode != 0 or not assembly.exists():
        raise RuntimeError(f"SPAdes failed for {output_dir.parent.name}")
    return assembly


def _map_to_draft(
    assembly: Path, r1: Path, r2: Path, bam: Path, *, threads: int, aligner: str,
    samtools: str,
) -> None:
    if bam.exists() and Path(str(bam) + ".bai").exists():
        return
    index = Path(str(assembly) + (".bwt.2bit.64" if aligner == "bwa-mem2" else ".bwt"))
    if not index.exists():
        subprocess.run([aligner, "index", str(assembly)], check=True)
    log = bam.with_suffix(".mapping.log")
    with log.open("w") as log_handle:
        mapping = subprocess.Popen(
            [aligner, "mem", "-t", str(threads), str(assembly), str(r1), str(r2)],
            stdout=subprocess.PIPE, stderr=log_handle,
        )
        assert mapping.stdout is not None
        sorting = subprocess.run(
            [samtools, "sort", "-@", str(threads), "-o", str(bam), "-"],
            stdin=mapping.stdout, stdout=log_handle, stderr=log_handle,
        )
        mapping.stdout.close()
        mapping_code = mapping.wait()
    if mapping_code != 0 or sorting.returncode != 0:
        bam.unlink(missing_ok=True)
        raise RuntimeError(f"read mapping failed for {bam.parent.name}")
    subprocess.run([samtools, "index", str(bam)], check=True)


def _discover(
    assembly: Path, baits: Path, bait_lengths: dict[str, int], sample_dir: Path,
) -> dict[str, tuple[list[str], list[tuple[int, int]]]]:
    database = sample_dir / "draft_db"
    if not Path(str(database) + ".nsq").exists():
        subprocess.run(
            ["makeblastdb", "-in", str(assembly), "-dbtype", "nucl", "-out", str(database)],
            check=True, capture_output=True, text=True,
        )
    output = sample_dir / "fixed_baits_vs_draft.tsv"
    fields = "qseqid sseqid pident length qstart qend sstart send bitscore"
    subprocess.run(
        [
            "blastn", "-query", str(baits), "-db", str(database), "-out", str(output),
            "-outfmt", f"6 {fields}", "-max_target_seqs", "20", "-perc_identity", "80",
        ], check=True,
    )
    candidates: dict[str, list[dict]] = {"23S": [], "gyrB": []}
    with output.open() as handle:
        for line in handle:
            qid, sid, pid, length, qs, qe, ss, se, bitscore = line.rstrip().split("\t")
            locus = qid.split("|", 1)[0]
            if locus not in candidates or int(length) < max(200, int(0.10 * bait_lengths[locus])):
                continue
            candidates[locus].append({
                "contig": sid,
                "subject_start": min(int(ss), int(se)),
                "subject_end": max(int(ss), int(se)),
                "query_start": min(int(qs), int(qe)),
                "query_end": max(int(qs), int(qe)),
                "bitscore": float(bitscore),
                "identity": float(pid),
            })
    discovered = {}
    for locus, hits in candidates.items():
        retained: list[dict] = []
        for hit in sorted(hits, key=lambda item: item["bitscore"], reverse=True):
            overlaps = any(
                old["contig"] == hit["contig"]
                and hit["subject_start"] <= old["subject_end"]
                and old["subject_start"] <= hit["subject_end"]
                for old in retained
            )
            if not overlaps:
                retained.append(hit)
        if not retained:
            raise RuntimeError(f"no qualifying {locus} hit in {assembly}")
        retained.sort(key=lambda item: (item["contig"], item["subject_start"]))
        regions = [
            f"{hit['contig']}:{hit['subject_start']}-{hit['subject_end']}" for hit in retained
        ]
        spans = [(hit["query_start"], hit["query_end"]) for hit in retained]
        discovered[locus] = (regions, spans)
    return discovered


def _require_legacy_depth_compatibility(baseline: Path, candidate: Path) -> None:
    """Fail if the graph extension changes any archived depth-result field."""

    require_files([baseline, candidate])
    with baseline.open(newline="") as handle:
        baseline_rows = list(csv.DictReader(handle, delimiter="\t"))
    with candidate.open(newline="") as handle:
        candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
    for label, rows in (("baseline", baseline_rows), ("candidate", candidate_rows)):
        if not rows:
            raise RuntimeError(f"{label} low-copy result contains no rows")
        missing = [field for field in LEGACY_RESULT_FIELDS if field not in rows[0]]
        if missing:
            raise RuntimeError(f"{label} result omits legacy fields: {', '.join(missing)}")
    key = lambda row: (row["strain"], row["locus"])
    expected = {key(row): row for row in baseline_rows}
    observed = {key(row): row for row in candidate_rows}
    if expected.keys() != observed.keys():
        raise RuntimeError("candidate low-copy rows differ from the archived sample/locus set")
    mismatches = [
        f"{strain}/{locus}:{field}={expected[(strain, locus)][field]!r}"
        f"->{observed[(strain, locus)][field]!r}"
        for strain, locus in sorted(expected)
        for field in LEGACY_RESULT_FIELDS
        if expected[(strain, locus)][field] != observed[(strain, locus)][field]
    ]
    if mismatches:
        raise RuntimeError(
            "graph/depth consensus changed archived depth fields: "
            + "; ".join(mismatches[:10])
        )


def main() -> int:
    repo_root = REPO_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path,
        default=repo_root / "validation" / "manifests" / "low_copy_23S_prjna816422.tsv",
    )
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-gb", type=int, default=8)
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260817)
    args = parser.parse_args()
    if min(args.threads, args.memory_gb, args.n_boot) < 1:
        parser.error("threads, memory and bootstrap count must be positive")

    require_files([args.manifest])
    with args.manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if len(manifest) != 5:
        raise ValueError(f"frozen low-copy manifest must contain five samples, found {len(manifest)}")
    if [int(row["selection_order"]) for row in manifest] != list(range(1, 6)):
        raise ValueError("low-copy manifest selection order is not deterministic")

    for tool in ("spades.py", "blastn", "makeblastdb", "samtools"):
        if not shutil.which(tool):
            raise FileNotFoundError(f"{tool} is required on PATH")
    aligner = "bwa-mem2" if shutil.which("bwa-mem2") else "bwa"
    if not shutil.which(aligner):
        raise FileNotFoundError("bwa-mem2 or bwa is required on PATH")

    args.work_dir.mkdir(parents=True, exist_ok=True)
    reads_dir = args.work_dir / "reads"
    truth_dir = args.work_dir / "closed_truth"
    baits, bait_lengths = _prepare_fixed_baits(args.work_dir)
    truth_audit: list[dict] = []
    read_provenance: list[dict] = []
    result_rows: list[dict] = []
    result_records: list[dict] = []
    provenance_inputs: list[Path] = [args.manifest, baits]

    for row in manifest:
        truth, gbff = _audit_truth(row, truth_dir)
        truth_audit.append(truth)
        provenance_inputs.append(gbff)
        r1, r2, read_records = _prepare_reads(row["illumina_run"], reads_dir)
        read_provenance.extend(read_records)
        provenance_inputs.extend([r1, r2])
        sample_dir = args.work_dir / "samples" / row["strain"]
        sample_dir.mkdir(parents=True, exist_ok=True)
        assembly = _assemble_draft(
            r1, r2, sample_dir / "spades", threads=args.threads,
            memory_gb=args.memory_gb, spades="spades.py",
        )
        gfa = sample_dir / "spades" / "assembly_graph_after_simplification.gfa"
        require_files([gfa])
        graph = parse_gfa(str(gfa))
        provenance_inputs.append(gfa)
        bam = sample_dir / "reads_vs_illumina_draft.sorted.bam"
        _map_to_draft(
            assembly, r1, r2, bam, threads=args.threads, aligner=aligner,
            samtools="samtools",
        )
        discovered = _discover(assembly, baits, bait_lengths, sample_dir)
        sample_record = {"strain": row["strain"], "illumina_run": row["illumina_run"]}
        for locus, truth_copies in (("23S", 2), ("gyrB", 1)):
            regions, spans = discovered[locus]
            result = estimate_locus_copy_number(
                bam,
                regions,
                gene_spans=spans,
                gene_length=bait_lengths[locus],
                max_backbone_bootstrap_blocks=512,
                n_boot=args.n_boot,
                seed=args.seed,
                samtools="samtools",
            )
            discovery = discover_graph_target_hits(
                graph,
                baits,
                bait_id=f"{locus}|source={FIXED_REFERENCE}",
                blastn="blastn",
            )
            graph_result = infer_graph_contexts(
                graph,
                discovery.hits,
                bait_length=discovery.bait_length,
            )
            consensus = reconcile_copy_number(result, graph_result)
            absolute_error = (
                abs(result.dosage_estimate - truth_copies)
                if result.dosage_estimate is not None else None
            )
            copy_number_absolute_error = (
                abs(consensus.copy_number_call - truth_copies)
                if consensus.copy_number_call is not None else None
            )
            row_result = {
                "strain": row["strain"],
                "biosample": row["biosample"],
                "illumina_run": row["illumina_run"],
                "assembly_accession": row["assembly_accession"],
                "locus": locus,
                "truth_copies": truth_copies,
                "assembly_region_count": result.assembly_region_count,
                "ratio": result.ratio,
                "ratio_ci_low": result.ratio_ci_low,
                "ratio_ci_high": result.ratio_ci_high,
                "depth_call": result.depth_call,
                "dosage_estimate": result.dosage_estimate,
                "dosage_status": result.dosage_status,
                "dosage_flags": ";".join(result.dosage_flags),
                "absolute_error": absolute_error,
                "ambiguity_index": result.ambiguity_index,
                "gene_total_aligned_bp": result.gene_total_aligned_bp,
                "gene_union_covered_bp": result.gene_union_covered_bp,
                "gene_overlap_bp": result.gene_overlap_bp,
                "locus_regions": ";".join(result.locus_regions),
                "region_gene_spans": ";".join(f"{a}-{b}" for a, b in result.region_gene_spans),
                "graph_target_segments": ";".join(graph_result.target_segments),
                "graph_bait_covered_bp": graph_result.bait_covered_bp,
                "graph_bait_coverage_fraction": graph_result.bait_coverage_fraction,
                "left_context_count": graph_result.left_context_count,
                "right_context_count": graph_result.right_context_count,
                "graph_context_count": graph_result.context_count,
                "graph_context_lower_bound": graph_result.context_lower_bound,
                "graph_context_status": graph_result.status,
                "graph_context_flags": ";".join(graph_result.flags),
                "copy_number_call": consensus.copy_number_call,
                "copy_number_kind": consensus.copy_number_kind,
                "copy_number_method": consensus.copy_number_method,
                "copy_number_status": consensus.copy_number_status,
                "copy_number_lower_bound": consensus.copy_number_lower_bound,
                "copy_number_upper_bound": consensus.copy_number_upper_bound,
                "consensus_flags": ";".join(consensus.consensus_flags),
                "copy_number_absolute_error": copy_number_absolute_error,
            }
            result_rows.append(row_result)
            sample_record[locus] = result.as_row()
            sample_record[locus].update({
                "graph_target_discovery": discovery.as_dict(),
                "graph_context": graph_result.as_dict(),
                "copy_number_consensus": consensus.as_dict(),
            })
        result_records.append(sample_record)
        print(
            f"{row['strain']}\t23S={sample_record['23S']['depth_call']}"
            f"/copy={sample_record['23S']['copy_number_consensus']['copy_number_call']}\t"
            f"gyrB={sample_record['gyrB']['depth_call']}"
            f"/copy={sample_record['gyrB']['copy_number_consensus']['copy_number_call']}",
            flush=True,
        )

    results_dir = args.work_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    write_json(results_dir / "closed_truth_audit.json", truth_audit)
    write_json(results_dir / "read_provenance.json", read_provenance)
    write_json(results_dir / "low_copy_results.json", result_records)
    result_table = results_dir / "low_copy_results.tsv"
    with result_table.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(result_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(result_rows)
    _require_legacy_depth_compatibility(
        repo_root / "validation" / "low-copy-23S" / "results.tsv",
        result_table,
    )

    rrna = [row for row in result_rows if row["locus"] == "23S"]
    controls = [row for row in result_rows if row["locus"] == "gyrB"]
    two_copy_rejections = sum(row["depth_call"] == "MULTICOPY_DEPTH" for row in rrna)
    control_false_positives = sum(row["depth_call"] == "MULTICOPY_DEPTH" for row in controls)
    exact_two_copy_consensus = sum(
        row["copy_number_call"] == 2
        and row["copy_number_kind"] == "INTEGER_CONTEXT_COUNT"
        and row["copy_number_status"] == "SUPPORTED"
        for row in rrna
    )
    exact_single_copy_consensus = sum(
        row["copy_number_call"] == 1
        and row["copy_number_kind"] == "INTEGER_CONTEXT_COUNT"
        and row["copy_number_status"] == "SUPPORTED"
        for row in controls
    )
    evaluation = {
        "schema_version": "1.0",
        "samples": len(manifest),
        "two_copy_rejections": two_copy_rejections,
        "two_copy_sensitivity": two_copy_rejections / len(rrna),
        "single_copy_control_false_positives": control_false_positives,
        "criterion_two_copy_sensitivity_ge_0_80": two_copy_rejections >= 4,
        "criterion_zero_control_false_positives": control_false_positives == 0,
        "criterion_legacy_depth_fields_unchanged": True,
        "exact_two_copy_consensus_calls": exact_two_copy_consensus,
        "exact_single_copy_control_consensus_calls": exact_single_copy_consensus,
        "criterion_all_two_copy_consensus_calls_correct": exact_two_copy_consensus == len(rrna),
        "criterion_all_single_copy_control_consensus_calls_correct": (
            exact_single_copy_consensus == len(controls)
        ),
    }
    evaluation["all_primary_criteria_passed"] = (
        evaluation["criterion_two_copy_sensitivity_ge_0_80"]
        and evaluation["criterion_zero_control_false_positives"]
        and evaluation["criterion_legacy_depth_fields_unchanged"]
        and evaluation["criterion_all_two_copy_consensus_calls_correct"]
        and evaluation["criterion_all_single_copy_control_consensus_calls_correct"]
    )
    write_json(results_dir / "evaluation.json", evaluation)

    provenance_inputs.extend([result_table, results_dir / "closed_truth_audit.json"])
    provenance = build_provenance(
        workflow="external-two-copy-23S",
        repo_root=repo_root,
        inputs=provenance_inputs,
        parameters={
            "selection": "first five eligible BioSamples by numeric strain suffix",
            "fixed_bait_reference": FIXED_REFERENCE,
            "draft_assembler": "SPAdes --isolate",
            "threads": args.threads,
            "memory_gb": args.memory_gb,
            "n_boot": args.n_boot,
            "max_backbone_bootstrap_blocks": 512,
            "seed": args.seed,
            "call_threshold": 1.5,
            "min_bq": 20,
            "min_mq": 20,
            "graph_min_identity": 80.0,
            "graph_min_query_coverage": 0.90,
            "graph_min_context_bp": "max(100, 2 * largest GFA overlap)",
            "graph_max_context_nodes": 8,
            "graph_max_context_paths": 64,
            "graph_max_context_bp": 5000,
        },
        tools={
            "python": [sys.executable, "--version"],
            "biopython": [sys.executable, "-c", "import Bio; print(Bio.__version__)"],
            "spades": ["spades.py", "--version"],
            "blastn": ["blastn", "-version"],
            "aligner": [aligner, "version"],
            "samtools": ["samtools", "--version"],
        },
    )
    provenance["complete"] = True
    provenance["primary_criteria_passed"] = evaluation["all_primary_criteria_passed"]
    write_json(results_dir / "workflow_provenance.json", provenance)
    return 0 if evaluation["all_primary_criteria_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
