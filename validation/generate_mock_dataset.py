#!/usr/bin/env python3
"""Generate a compact deterministic truth set for Locus-Recon validation."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import random
import shutil
from pathlib import Path
from typing import Iterable, List, Tuple

SEED = 731_921
READ_LENGTH = 150
INSERT_SIZE = 350
ERROR_RATE = 0.001
LOCUS_NAME = "mockLocus"
STOP_CODONS = {"TAA", "TAG", "TGA"}
DNA = "ACGT"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(DNA) for _ in range(length))


def build_allele(rng: random.Random, codons: int = 169) -> str:
    allowed = [
        a + b + c
        for a in DNA for b in DNA for c in DNA
        if a + b + c not in STOP_CODONS
    ]
    sequence = ["ATG"]
    sequence.extend(rng.choice(allowed) for _ in range(codons - 1))
    return "".join(sequence)


def mutate_codons(
    sequence: str,
    codon_positions: Iterable[int],
    rng: random.Random,
) -> Tuple[str, List[int]]:
    codons = [sequence[index:index + 3] for index in range(0, len(sequence), 3)]
    changed_bases = []
    for codon_position in codon_positions:
        original = codons[codon_position]
        alternatives = [
            a + b + c
            for a in DNA for b in DNA for c in DNA
            if a + b + c not in STOP_CODONS and a + b + c != original
        ]
        replacement = rng.choice(alternatives)
        codons[codon_position] = replacement
        changed_bases.extend(
            codon_position * 3 + offset + 1
            for offset, (old, new) in enumerate(zip(original, replacement))
            if old != new
        )
    return "".join(codons), changed_bases


def introduce_errors(sequence: str, rng: random.Random) -> str:
    bases = list(sequence)
    for index, base in enumerate(bases):
        if rng.random() < ERROR_RATE:
            bases[index] = rng.choice([candidate for candidate in DNA if candidate != base])
    return "".join(bases)


def simulate_pairs(
    genome: str,
    coverage: float,
    seed: int,
    prefix: str,
) -> List[Tuple[str, str, str]]:
    """Simulate evenly distributed, bidirectionally oriented paired reads."""
    rng = random.Random(seed)
    number_of_pairs = max(1, math.ceil(coverage * len(genome) / (2 * READ_LENGTH)))
    maximum_start = len(genome) - INSERT_SIZE
    pairs = []
    for pair_index in range(number_of_pairs):
        if number_of_pairs == 1:
            start = maximum_start // 2
        else:
            expected = pair_index * maximum_start / (number_of_pairs - 1)
            start = int(round(expected + rng.randint(-8, 8)))
            start = max(0, min(maximum_start, start))
        fragment = genome[start:start + INSERT_SIZE]
        if pair_index % 2:
            read_1 = reverse_complement(fragment[-READ_LENGTH:])
            read_2 = fragment[:READ_LENGTH]
        else:
            read_1 = fragment[:READ_LENGTH]
            read_2 = reverse_complement(fragment[-READ_LENGTH:])
        name = f"{prefix}_{pair_index:06d}_{start}"
        pairs.append((
            name,
            introduce_errors(read_1, rng),
            introduce_errors(read_2, rng),
        ))
    return pairs


def simulate_mixture_pairs(
    genome_a: str,
    genome_b: str,
    total_coverage: float,
    fraction_b: float,
    seed: int,
    prefix: str,
) -> List[Tuple[str, str, str]]:
    pairs_a = simulate_pairs(
        genome_a, total_coverage * (1.0 - fraction_b), seed + 11, prefix + "_A"
    )
    pairs_b = simulate_pairs(
        genome_b, total_coverage * fraction_b, seed + 29, prefix + "_B"
    ) if fraction_b else []
    pairs = pairs_a + pairs_b
    random.Random(seed + 47).shuffle(pairs)
    return pairs


def write_fasta(path: Path, records: Iterable[Tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def _gzip_text_writer(path: Path):
    raw_handle = path.open("wb")
    gzip_handle = gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0)
    return raw_handle, gzip_handle, io.TextIOWrapper(gzip_handle, encoding="ascii", newline="\n")


def write_fastq_pairs(
    r1_path: Path,
    r2_path: Path,
    pairs: Iterable[Tuple[str, str, str]],
) -> None:
    raw_1, gzip_1, handle_1 = _gzip_text_writer(r1_path)
    raw_2, gzip_2, handle_2 = _gzip_text_writer(r2_path)
    try:
        for name, read_1, read_2 in pairs:
            handle_1.write(f"@{name}/1\n{read_1}\n+\n{'I' * len(read_1)}\n")
            handle_2.write(f"@{name}/2\n{read_2}\n+\n{'I' * len(read_2)}\n")
    finally:
        handle_1.close()
        handle_2.close()
        if not gzip_1.closed:
            gzip_1.close()
        if not gzip_2.closed:
            gzip_2.close()
        if not raw_1.closed:
            raw_1.close()
        if not raw_2.closed:
            raw_2.close()


def write_long_reads(path: Path, genome: str, locus_start: int) -> None:
    rng = random.Random(SEED + 991)
    raw, gzip_handle, handle = _gzip_text_writer(path)
    try:
        for index in range(12):
            start = max(0, locus_start - 700 + index * 30)
            sequence = genome[start:start + 1900]
            if index % 2:
                sequence = reverse_complement(sequence)
            sequence = introduce_errors(sequence, rng)
            handle.write(
                f"@mock_long_{index:03d}\n{sequence}\n+\n{'I' * len(sequence)}\n"
            )
    finally:
        handle.close()
        if not gzip_handle.closed:
            gzip_handle.close()
        if not raw.closed:
            raw.close()


def sha256(path: Path) -> str:
    """Return the SHA-256 checksum of one generated benchmark file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fragmented_assembly(
    genome: str,
    locus_start: int,
    locus_length: int,
    split_fraction: float,
) -> List[Tuple[str, str]]:
    breakpoint = locus_start + int(round(locus_length * split_fraction))
    return [
        ("fragment_left", genome[:breakpoint]),
        ("fragment_right", genome[breakpoint:]),
    ]


def build_dataset(output_dir: Path, force: bool = False) -> None:
    if output_dir.exists():
        if not force:
            raise FileExistsError(
                f"{output_dir} already exists; use --force to regenerate it."
            )
        shutil.rmtree(output_dir)

    assemblies_dir = output_dir / "assemblies"
    reads_dir = output_dir / "reads"
    truth_dir = output_dir / "truth"
    for directory in (assemblies_dir, reads_dir, truth_dir):
        directory.mkdir(parents=True, exist_ok=True)

    rng = random.Random(SEED)
    allele_a = build_allele(rng)
    allele_b, allele_b_changes = mutate_codons(
        allele_a, [7, 18, 31, 44, 58, 73, 89, 104, 121, 139, 153, 164], rng
    )
    allele_c, allele_c_changes = mutate_codons(
        allele_a, [4, 13, 25, 38, 51, 65, 80, 96, 112, 128, 145, 160], rng
    )
    paralog, paralog_changes = mutate_codons(
        allele_a, [20, 55, 90, 125, 155], rng
    )

    left = random_dna(1250, rng)
    right = random_dna(1250, rng)
    spacer = random_dna(850, rng)
    genome_a = left + allele_a + right
    genome_b = left + allele_b + right
    paralog_genome = left + allele_a + spacer + paralog + right
    negative_genome = random_dna(len(genome_a), random.Random(SEED + 8_000))
    locus_start = len(left)

    write_fasta(output_dir / "mock_locus_alleles.fasta", [
        (f"{LOCUS_NAME}_A", allele_a),
        (f"{LOCUS_NAME}_B", allele_b),
        (f"{LOCUS_NAME}_C", allele_c),
    ])
    write_fasta(truth_dir / "truth_alleles.fasta", [
        (f"{LOCUS_NAME}_A", allele_a),
        (f"{LOCUS_NAME}_B", allele_b),
        (f"{LOCUS_NAME}_paralog", paralog),
    ])
    write_fasta(truth_dir / "complete_genome_A.fasta", [("complete_A", genome_a)])
    write_fasta(truth_dir / "complete_genome_B.fasta", [("complete_B", genome_b)])
    write_fasta(
        truth_dir / "complete_genome_paralog.fasta",
        [("complete_A_plus_paralog", paralog_genome)],
    )
    write_fasta(
        truth_dir / "complete_genome_negative.fasta",
        [("complete_negative", negative_genome)],
    )
    write_long_reads(truth_dir / "mock_long_reads_A.fastq.gz", genome_a, locus_start)

    write_fasta(assemblies_dir / "intact_A.fasta", [("chromosome", genome_a)])
    write_fasta(
        assemblies_dir / "fragmented_A.fasta",
        fragmented_assembly(genome_a, locus_start, len(allele_a), 0.50),
    )
    write_fasta(
        assemblies_dir / "fragmented_A_terminal.fasta",
        fragmented_assembly(genome_a, locus_start, len(allele_a), 0.12),
    )
    target_start = locus_start
    paralog_start = len(left) + len(allele_a) + len(spacer)
    flank = 550
    write_fasta(assemblies_dir / "paralog_draft.fasta", [
        ("target_copy", paralog_genome[target_start - flank:target_start + len(allele_a) + flank]),
        (
            "paralog_copy",
            paralog_genome[
                paralog_start - flank:paralog_start + len(paralog) + flank
            ],
        ),
    ])
    write_fasta(
        assemblies_dir / "negative.fasta",
        [("negative_contig", negative_genome)],
    )

    samples = [
        {
            "sample_id": "baseline_intact_A", "category": "baseline",
            "assembly": "assemblies/intact_A.fasta", "genome": genome_a,
            "depth": 40.0, "mixture_fraction": 0.0, "truth_allele": f"{LOCUS_NAME}_A",
            "secondary_allele": "",
            "expected": "exact reconstruction",
        },
        {
            "sample_id": "fragmented_mid_A", "category": "fragmented",
            "assembly": "assemblies/fragmented_A.fasta", "genome": genome_a,
            "depth": 40.0, "mixture_fraction": 0.0, "truth_allele": f"{LOCUS_NAME}_A",
            "secondary_allele": "",
            "expected": "exact reconstruction",
        },
        {
            "sample_id": "fragmented_terminal_A", "category": "fragmented",
            "assembly": "assemblies/fragmented_A_terminal.fasta", "genome": genome_a,
            "depth": 40.0, "mixture_fraction": 0.0, "truth_allele": f"{LOCUS_NAME}_A",
            "secondary_allele": "",
            "expected": "exact reconstruction",
        },
    ]
    for depth in (3.0, 5.0, 10.0, 20.0):
        samples.append({
            "sample_id": f"depth_{int(depth):02d}x_A", "category": "depth_series",
            "assembly": "assemblies/fragmented_A.fasta", "genome": genome_a,
            "depth": depth, "mixture_fraction": 0.0,
            "truth_allele": f"{LOCUS_NAME}_A", "expected": "sensitivity series",
            "secondary_allele": "",
        })
    for fraction in (0.05, 0.10, 0.20, 0.50):
        samples.append({
            "sample_id": f"mixture_{int(fraction * 100):02d}pct_B",
            "category": "mixture",
            "assembly": "assemblies/fragmented_A.fasta", "genome": genome_a,
            "genome_b": genome_b, "depth": 60.0, "mixture_fraction": fraction,
            "truth_allele": f"{LOCUS_NAME}_A", "expected": "mixture detected",
            "secondary_allele": f"{LOCUS_NAME}_B",
        })
    samples.extend([
        {
            "sample_id": "paralog_control", "category": "paralog",
            "assembly": "assemblies/paralog_draft.fasta", "genome": paralog_genome,
            "depth": 60.0, "mixture_fraction": 0.0,
            "truth_allele": f"{LOCUS_NAME}_A", "expected": "ambiguity flagged",
            "secondary_allele": "",
        },
        {
            "sample_id": "negative_control", "category": "negative",
            "assembly": "assemblies/negative.fasta", "genome": negative_genome,
            "depth": 40.0, "mixture_fraction": 0.0,
            "truth_allele": "", "expected": "no reconstruction",
            "secondary_allele": "",
        },
    ])

    for index, sample in enumerate(samples):
        sample_id = sample["sample_id"]
        if sample["category"] == "mixture":
            pairs = simulate_mixture_pairs(
                sample["genome"], sample["genome_b"], sample["depth"],
                sample["mixture_fraction"], SEED + 1000 + index * 101, sample_id,
            )
        else:
            pairs = simulate_pairs(
                sample["genome"], sample["depth"],
                SEED + 1000 + index * 101, sample_id,
            )
        r1_relative = f"reads/{sample_id}_R1.fastq.gz"
        r2_relative = f"reads/{sample_id}_R2.fastq.gz"
        write_fastq_pairs(output_dir / r1_relative, output_dir / r2_relative, pairs)
        sample["r1"] = r1_relative
        sample["r2"] = r2_relative
        sample["read_pairs"] = len(pairs)

    with (output_dir / "samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "assembly_path", "r1_path", "r2_path"])
        for sample in samples:
            writer.writerow([
                sample["sample_id"], sample["assembly"], sample["r1"], sample["r2"]
            ])

    with (truth_dir / "expected_results.tsv").open("w", newline="") as handle:
        columns = [
            "sample_id", "category", "truth_allele", "secondary_allele", "nominal_depth",
            "mixture_fraction", "expected_outcome", "read_pairs",
        ]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for sample in samples:
            writer.writerow({
                "sample_id": sample["sample_id"],
                "category": sample["category"],
                "truth_allele": sample["truth_allele"],
                "secondary_allele": sample["secondary_allele"],
                "nominal_depth": sample["depth"],
                "mixture_fraction": sample["mixture_fraction"],
                "expected_outcome": sample["expected"],
                "read_pairs": sample["read_pairs"],
            })

    manifest_samples = [
        {key: value for key, value in sample.items() if key not in {"genome", "genome_b"}}
        for sample in samples
    ]
    central_break = locus_start + int(round(len(allele_a) * 0.50))
    terminal_break = locus_start + int(round(len(allele_a) * 0.12))
    generated_hashes = {
        str(path.relative_to(output_dir)): sha256(path)
        for path in sorted(output_dir.rglob("*"))
        if path.is_file()
    }
    manifest = {
        "schema_version": "1.0",
        "generator": "validation/generate_mock_dataset.py",
        "seed": SEED,
        "locus": LOCUS_NAME,
        "allele_length": len(allele_a),
        "complete_genome_A_length": len(genome_a),
        "locus_start_0_based": locus_start,
        "locus_end_0_based_exclusive": locus_start + len(allele_a),
        "fragmentation_breaks_0_based": {
            "central": central_break,
            "near_terminal": terminal_break,
        },
        "mock_long_read_count": 12,
        "read_length": READ_LENGTH,
        "insert_size": INSERT_SIZE,
        "substitution_error_rate": ERROR_RATE,
        "allele_b_changed_positions_1_based": allele_b_changes,
        "allele_c_changed_positions_1_based": allele_c_changes,
        "paralog_changed_positions_1_based": paralog_changes,
        "truth_model": (
            "Complete synthetic chromosomes plus locus-spanning mock long reads "
            "act as the hybrid/long-read ground truth."
        ),
        "file_sha256_excluding_manifest": generated_hashes,
        "samples": manifest_samples,
    }
    with (output_dir / "dataset_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "mock_dataset",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    build_dataset(args.output_dir.resolve(), force=args.force)
    print(f"Mock validation dataset written to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
