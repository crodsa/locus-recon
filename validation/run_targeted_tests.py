#!/usr/bin/env python3
"""Targeted tests of Locus-Recon with simulated reads (protocol section 3).

This runner lives in ``validation/`` beside ``generate_mock_dataset.py`` and
``common.py`` and imports both.  It never edits the package: it calls the
installed command-line programs (``locus-recon``, ``locus-recon-depth-ratio``,
``locus-recon-copy-number``) with default options and, for the block-size
series only, the library function ``estimate_locus_copy_number``.

Tests (``validation/external-evaluation/PROTOCOL.md``, section 3):

1. mixtures of two alleles that differ at few sites, and pure controls;
2. strand-bias alerts in pure samples by locus length and depth;
3. terminal substitutions, indels and pseudogenes;
4. graph contexts with known copy number;
5. overlapping discovery spans with one and two copies;
   plus bootstrap block sizes 100, 250, 500 and 1,000 bp for tests 4 and 5.

Every simulated sample is: one or two synthetic genomes -> paired reads with
substitution errors -> SPAdes ``--isolate`` draft -> Locus-Recon.  All seeds are
derived from ``BASE_SEED`` and the run identifier, so a run is reproduced by
its identifier alone.

Example::

    python validation/run_targeted_tests.py --work-dir /scratch/targeted \
        --out-dir validation/targeted-tests --jobs 8 --threads-per-run 1
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time
import traceback
import zlib
from collections import Counter, OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

VALIDATION_DIR = Path(__file__).resolve().parent
if str(VALIDATION_DIR) not in sys.path:
    sys.path.insert(0, str(VALIDATION_DIR))

import numpy as np

from common import command_version, sha256_file, source_sha256, write_json
from generate_mock_dataset import (
    ERROR_RATE,
    INSERT_SIZE,
    READ_LENGTH,
    STOP_CODONS,
    reverse_complement,
    write_fasta,
)

# --------------------------------------------------------------------------- #
# Design constants.  Fixed before any result was seen (see README.md).
# --------------------------------------------------------------------------- #
BASE_SEED = 20261118
EXPECTED_PACKAGE_DIGEST = (
    "2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0"
)
EXPECTED_VERSION = "1.0.1"

INSERT_SD = 30            # the generator uses a fixed 350-bp fragment
INSERT_MIN, INSERT_MAX = 200, 500
BACKGROUND_BP = 100_000   # random background, as run_constructed_copy_series.py
# Background arms of the copy-number tests.  "plain": the locus is the only
# repeat of a linear genome (the design written before the pilot).  "decoy": an
# unrelated 2,000-bp repeat is placed at the midpoint of the last two background
# segments, so that it interleaves with the locus copies and the assembly graph
# has more than one traversal, as the graph of a bacterial genome has (added
# after the pilot; see DEVIATIONS.md).
DECOY_REPEAT_BP = 2000
T4_BACKGROUNDS = ("plain", "decoy")
T5_BACKGROUND = "decoy"
LOCUS = "tlocus"
BAIT_ID = "tlocus_ref"

# Test 1
T1_LOCUS_BP = 1500
T1_SITES = (1, 2, 5, 10)
T1_FRACTIONS = (5, 10, 20, 35, 50)
T1_DEPTH = 60
T1_REPLICATES = 10
T1_PURE_CONTROLS = 20
# Test 2 (500 and 7,000 bp are not multiples of three; the coding alleles are
# 501 and 6,999 bp)
T2_LENGTHS = (501, 1500, 3000, 6999)
T2_NOMINAL = {501: 500, 1500: 1500, 3000: 3000, 6999: 7000}
T2_DEPTHS = (20, 60, 200, 1000)
T2_REPLICATES = 10
# Test 3
T3_LOCUS_BP = 1500
T3_SCENARIOS = (
    "control_same_as_bait", "sub_ends_5bp", "sub_ends_20bp", "sub_ends_50bp",
    "del_3bp", "ins_12bp", "del_1bp", "compensating_indels_30bp",
    "premature_stop",
)
T3_DEPTHS = (10, 20, 60)
T3_REPLICATES = 5
T3_INTERIOR_MARGIN = 200  # indels and stops are placed at least this far from an end
# Test 4
T4_LOCUS_BP = 3000
T4_DEPTH = 60
T4_REPLICATES = 5
T4_DIVERGENCE = 0.01
T4_FLANK_SITE_DISTANCE = 300   # strain-mixture site: bp upstream of the locus
T4_FLANK_REPEAT_BP = 1000
T4_FLANK_REPEAT_SPACER = 100   # unique bp between the flank repeat and the locus
T4_SCENARIOS = (
    # name, copies, layout, divergent
    ("c1_control", 1, "dispersed", False),
    ("c2_dispersed_identical", 2, "dispersed", False),
    ("c2_dispersed_div1pct", 2, "dispersed", True),
    ("c3_dispersed_identical", 3, "dispersed", False),
    ("c3_dispersed_div1pct", 3, "dispersed", True),
    ("c4_dispersed_identical", 4, "dispersed", False),
    ("c4_dispersed_div1pct", 4, "dispersed", True),
    ("c2_sharedflank500_identical", 2, "sharedflank500", False),
    ("c2_sharedflank500_div1pct", 2, "sharedflank500", True),
    ("c2_sharedflank2000_identical", 2, "sharedflank2000", False),
    ("c2_sharedflank2000_div1pct", 2, "sharedflank2000", True),
    ("c2_tandem_identical", 2, "tandem", False),
    ("c2_tandem_div1pct", 2, "tandem", True),
    ("c1_flankrepeat1kb", 1, "flankrepeat", False),
    ("c1_strainmix_flanksite", 1, "strainmix", False),
    ("c2_strainmix_flanksite_identical", 2, "strainmix", False),
    ("c2_strainmix_flanksite_div1pct", 2, "strainmix", True),
)
# Test 5
T5_LOCUS_BP = 3000
T5_DEPTH = 60
T5_REPLICATES = 5
T5_OVERLAPS = (0, 10, 25, 45, 65, 83, 90)
T5_COPIES = (1, 2)
# Block-size series (tests 4 and 5)
BLOCK_SIZES = (100, 250, 500, 1000)
BOOTSTRAP_SEED = 20260817   # package default
BOOTSTRAP_RESAMPLES = 1000  # package default

DNA = "ACGT"
_DNA_BYTES = np.frombuffer(b"ACGT", dtype="S1")
_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
SENSE_CODONS = [
    a + b + c for a in DNA for b in DNA for c in DNA if a + b + c not in STOP_CODONS
]
FLAG_NAME = re.compile(r"(?:^|; )([A-Z][A-Z0-9_]+)(?= \(|;|$)")


# --------------------------------------------------------------------------- #
# Seeds and sequences
# --------------------------------------------------------------------------- #
def seed_for(run_id: str, purpose: str) -> int:
    """Deterministic 63-bit seed from the base seed, the run and a purpose."""
    digest = hashlib.sha256(f"{BASE_SEED}|{run_id}|{purpose}".encode()).digest()
    return int.from_bytes(digest[:8], "big") >> 1


def random_dna(length: int, rng: np.random.Generator) -> str:
    return _DNA_BYTES[rng.integers(0, 4, size=length)].tobytes().decode("ascii")


def build_allele(length: int, rng: np.random.Generator) -> str:
    """ATG followed by sense codons, as ``generate_mock_dataset.build_allele``."""
    if length % 3:
        raise ValueError("coding alleles must be a multiple of three")
    picks = rng.integers(0, len(SENSE_CODONS), size=length // 3 - 1)
    return "ATG" + "".join(SENSE_CODONS[i] for i in picks)


def substitute(sequence: str, positions, rng: np.random.Generator):
    """Substitute 0-based positions without creating a stop codon in frame +1.

    Returns the new sequence and ``[(position_1_based, ref, alt), ...]``.
    """
    bases = list(sequence)
    changes = []
    for position in sorted(int(p) for p in positions):
        original = bases[position]
        alternatives = [b for b in DNA if b != original]
        order = rng.permutation(3)
        chosen = None
        codon_start = position - position % 3
        for index in order:
            trial = bases[codon_start:codon_start + 3]
            trial[position - codon_start] = alternatives[index]
            if len(trial) < 3 or "".join(trial) not in STOP_CODONS:
                chosen = alternatives[index]
                break
        if chosen is None:  # cannot happen for a sense codon; kept explicit
            raise RuntimeError(f"no stop-free substitution at {position}")
        bases[position] = chosen
        changes.append((position + 1, original, chosen))
    return "".join(bases), changes


def frame_facts(sequence: str) -> dict:
    """Stops in frame +1 and the six-frame minimum, with the package's counter."""
    from locus_recon.qc import _count_internal_stops, count_internal_stops_six_frames
    minimum, label, _counts = count_internal_stops_six_frames(sequence)
    return {
        "length_mod3": len(sequence) % 3,
        "internal_stops_frame_plus1": _count_internal_stops(sequence, 0),
        "internal_stops_min_six_frames": minimum,
        "min_stop_frame": label,
    }


# --------------------------------------------------------------------------- #
# Read simulation: the error model of generate_mock_dataset.py (substitutions
# at ERROR_RATE, quality 'I'), with uniformly random fragment starts as in
# workflows/run_constructed_copy_series.py, random orientation and a normal
# fragment length (mean INSERT_SIZE, sd INSERT_SD).
# --------------------------------------------------------------------------- #
def _encode(genome: str) -> np.ndarray:
    table = np.full(256, 255, dtype=np.uint8)
    for base, code in _CODE.items():
        table[ord(base)] = code
    encoded = table[np.frombuffer(genome.encode("ascii"), dtype=np.uint8)]
    if (encoded == 255).any():
        raise ValueError("genome holds a symbol other than A, C, G, T")
    return encoded


def simulate_pairs(genome: str, coverage: float, rng: np.random.Generator):
    """Return (R1 codes, R2 codes) arrays of shape (pairs, READ_LENGTH)."""
    codes = _encode(genome)
    pairs = max(1, math.ceil(coverage * len(codes) / (2 * READ_LENGTH)))
    inserts = np.clip(
        np.rint(rng.normal(INSERT_SIZE, INSERT_SD, size=pairs)), INSERT_MIN, INSERT_MAX
    ).astype(np.int64)
    starts = np.floor(rng.random(pairs) * (len(codes) - inserts + 1)).astype(np.int64)
    offsets = np.arange(READ_LENGTH, dtype=np.int64)
    left = codes[starts[:, None] + offsets]
    right = 3 - codes[(starts + inserts - 1)[:, None] - offsets]  # reverse complement
    flip = rng.random(pairs) < 0.5
    read_1 = np.where(flip[:, None], right, left)
    read_2 = np.where(flip[:, None], left, right)
    for reads in (read_1, read_2):
        errors = rng.random(reads.shape) < ERROR_RATE
        shift = rng.integers(1, 4, size=int(errors.sum()), dtype=np.uint8)
        reads[errors] = (reads[errors] + shift) % 4
    return read_1.astype(np.uint8), read_2.astype(np.uint8)


def simulate_sample(components, depth: float, rng: np.random.Generator):
    """Reads of one sample: ``components`` is ``[(genome, fraction), ...]``.

    As ``simulate_mixture_pairs`` of the generator: each component is sampled at
    ``depth * fraction`` and the pairs are shuffled together.
    """
    parts_1, parts_2 = [], []
    for genome, fraction in components:
        if fraction <= 0:
            continue
        read_1, read_2 = simulate_pairs(genome, depth * fraction, rng)
        parts_1.append(read_1)
        parts_2.append(read_2)
    read_1 = np.concatenate(parts_1)
    read_2 = np.concatenate(parts_2)
    order = rng.permutation(len(read_1))
    return read_1[order], read_2[order]


def write_fastq_pair(r1_path: Path, r2_path: Path, read_1, read_2, prefix: str) -> None:
    quality = "I" * READ_LENGTH
    lookup = np.frombuffer(b"ACGT", dtype=np.uint8)
    for path, reads, mate in ((r1_path, read_1, 1), (r2_path, read_2, 2)):
        text = lookup[reads].tobytes().decode("ascii")
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0,
                               compresslevel=6) as gz:
                buffer = io.StringIO()
                for index in range(len(reads)):
                    sequence = text[index * READ_LENGTH:(index + 1) * READ_LENGTH]
                    buffer.write(f"@{prefix}_{index:07d}/{mate}\n{sequence}\n+\n{quality}\n")
                    if index % 20000 == 19999:
                        gz.write(buffer.getvalue().encode("ascii"))
                        buffer = io.StringIO()
                gz.write(buffer.getvalue().encode("ascii"))


# --------------------------------------------------------------------------- #
# Manifest of runs
# --------------------------------------------------------------------------- #
def build_manifest(replicate_offset: int = 0, replicate_cap: int | None = None) -> list[dict]:
    """Every run of the five tests.  ``replicate_cap`` is for pilot runs only."""
    def reps(n):
        count = n if replicate_cap is None else min(n, replicate_cap)
        return [replicate_offset + i for i in range(1, count + 1)]

    runs: list[dict] = []
    for sites in T1_SITES:
        for fraction in T1_FRACTIONS:
            for rep in reps(T1_REPLICATES):
                runs.append({
                    "run_id": f"t1_s{sites:02d}_f{fraction:02d}_r{rep:02d}", "test": 1,
                    "scenario": "mixture", "sites": sites, "minor_fraction_pct": fraction,
                    "depth": T1_DEPTH, "locus_bp": T1_LOCUS_BP, "replicate": rep,
                })
    for rep in reps(T1_PURE_CONTROLS):
        runs.append({
            "run_id": f"t1_pure_r{rep:02d}", "test": 1, "scenario": "pure_control",
            "sites": 0, "minor_fraction_pct": 0, "depth": T1_DEPTH,
            "locus_bp": T1_LOCUS_BP, "replicate": rep,
        })
    for length in T2_LENGTHS:
        for depth in T2_DEPTHS:
            for rep in reps(T2_REPLICATES):
                runs.append({
                    "run_id": f"t2_L{T2_NOMINAL[length]:04d}_d{depth:04d}_r{rep:02d}",
                    "test": 2, "scenario": "pure", "depth": depth, "locus_bp": length,
                    "nominal_locus_bp": T2_NOMINAL[length], "replicate": rep,
                })
    for scenario in T3_SCENARIOS:
        for depth in T3_DEPTHS:
            for rep in reps(T3_REPLICATES):
                runs.append({
                    "run_id": f"t3_{scenario}_d{depth:02d}_r{rep:02d}", "test": 3,
                    "scenario": scenario, "depth": depth, "locus_bp": T3_LOCUS_BP,
                    "replicate": rep,
                })
    for background in T4_BACKGROUNDS:
        prefix = "t4" if background == "plain" else "t4d"
        for scenario, copies, layout, divergent in T4_SCENARIOS:
            for rep in reps(T4_REPLICATES):
                runs.append({
                    "run_id": f"{prefix}_{scenario}_r{rep:02d}", "test": 4,
                    "scenario": scenario, "background": background,
                    "true_copies": copies, "layout": layout, "divergent": divergent,
                    "depth": T4_DEPTH, "locus_bp": T4_LOCUS_BP, "replicate": rep,
                })
    for copies in T5_COPIES:
        for overlap in T5_OVERLAPS:
            for rep in reps(T5_REPLICATES):
                runs.append({
                    "run_id": f"t5_c{copies}_ov{overlap:02d}_r{rep:02d}", "test": 5,
                    "scenario": f"c{copies}_overlap{overlap:02d}",
                    "arm": f"c{copies}_identical_cut",
                    "background": T5_BACKGROUND, "true_copies": copies, "divergent": False,
                    "overlap_pct": overlap, "depth": T5_DEPTH, "locus_bp": T5_LOCUS_BP,
                    "replicate": rep,
                })
    # Added arm (see DEVIATIONS.md): two copies that differ at 1% of sites, each
    # kept in part by the draft, so that the reads of both map uniquely.
    for overlap in T5_OVERLAPS:
        for rep in reps(T5_REPLICATES):
            runs.append({
                "run_id": f"t5_c2div_ov{overlap:02d}_r{rep:02d}", "test": 5,
                "scenario": f"c2div_overlap{overlap:02d}", "arm": "c2_divergent_partial",
                "background": "plain", "true_copies": 2, "divergent": True,
                "overlap_pct": overlap, "depth": T5_DEPTH, "locus_bp": T5_LOCUS_BP,
                "replicate": rep,
            })
    return runs


# --------------------------------------------------------------------------- #
# Truth construction
# --------------------------------------------------------------------------- #
def _segments(count: int, rng, run: dict | None = None) -> list[str]:
    """Background segments; in the decoy arm the last two carry the decoy repeat."""
    size = BACKGROUND_BP // count
    parts = [random_dna(size, rng) for _ in range(count)]
    if run is not None and run.get("background") == "decoy":
        decoy = random_dna(
            DECOY_REPEAT_BP, np.random.default_rng(seed_for(run["run_id"], "decoy")))
        for index in (len(parts) - 2, len(parts) - 1):
            middle = len(parts[index]) // 2
            parts[index] = parts[index][:middle] + decoy + parts[index][middle:]
    return parts


def _background_note(run: dict) -> str:
    if run.get("background") == "decoy":
        return (f"; background with an unrelated {DECOY_REPEAT_BP}-bp repeat in two copies "
                "that interleave with the locus")
    return "; background without any other repeat"


def build_case(run: dict) -> dict:
    """Genomes, bait and truth of one run.

    Returns ``components`` ``[(name, genome, fraction)]``, ``bait``, ``alleles``
    ``[(name, sequence)]`` (every true copy/allele of the locus in the sample)
    and a ``truth`` dict that is written to the manifest.
    """
    rng = np.random.default_rng(seed_for(run["run_id"], "genome"))
    length = run["locus_bp"]
    test = run["test"]
    truth: dict = {}

    if test in (1, 2, 3):
        left, right = _segments(2, rng)
        reference = build_allele(length, rng)
        bait = reference
        if test == 1 and run["scenario"] == "mixture":
            positions = rng.choice(length, size=run["sites"], replace=False)
            minor, changes = substitute(reference, positions, rng)
            fraction = run["minor_fraction_pct"] / 100.0
            components = [
                ("strain_major", left + reference + right, 1.0 - fraction),
                ("strain_minor", left + minor + right, fraction),
            ]
            alleles = [("major", reference), ("minor", minor)]
            truth.update({
                "truth": f"two strains; alleles differ at {run['sites']} site(s); "
                         f"minor fraction {run['minor_fraction_pct']}%; bait = major allele",
                "variant_sites_1based": ";".join(f"{p}{r}>{a}" for p, r, a in changes),
                "true_copies_per_genome": 1,
            })
        elif test in (1, 2):
            components = [("strain", left + reference + right, 1.0)]
            alleles = [("allele", reference)]
            truth.update({
                "truth": "pure sample; allele identical to the bait",
                "variant_sites_1based": "", "true_copies_per_genome": 1,
            })
        else:
            sample, description, events = _test3_allele(run["scenario"], reference, rng)
            components = [("strain", left + sample + right, 1.0)]
            alleles = [("allele", sample)]
            facts = frame_facts(sample)
            truth.update({
                "truth": description,
                "variant_sites_1based": events,
                "true_copies_per_genome": 1,
                "truth_length": len(sample),
                "truth_length_mod3": facts["length_mod3"],
                "bait_length_mod3": len(bait) % 3,
                "truth_internal_stops_frame_plus1": facts["internal_stops_frame_plus1"],
                "truth_internal_stops_min_six_frames": facts["internal_stops_min_six_frames"],
                "truth_frame_intact": int(
                    facts["length_mod3"] == len(bait) % 3
                    and facts["internal_stops_frame_plus1"] == 0
                ),
            })
        truth["locus_start_0based"] = len(left)
        return {"components": components, "bait": bait, "alleles": alleles, "truth": truth}

    if test == 5:
        copies = run["true_copies"]
        reference = build_allele(length, rng)
        parts = _segments(copies + 1, rng, run)
        copy_sequences = [reference] * copies
        changes = []
        if run.get("divergent"):
            positions = rng.choice(length, size=round(T4_DIVERGENCE * length), replace=False)
            other, changes = substitute(reference, positions, rng)
            copy_sequences = [reference, other]
        genome = parts[0]
        for index in range(copies):
            genome += copy_sequences[index] + parts[index + 1]
        if run.get("divergent"):
            description = (
                "2 copies with different flanks that differ at 1% of sites; the draft "
                "keeps copy 1 from the start of the query and copy 2 to its end, and the "
                f"two kept parts share {run['overlap_pct']}% of the query")
        else:
            description = (
                f"{copies} identical cop{'y' if copies == 1 else 'ies'} with different "
                "flanks; the draft contig that holds the locus is cut in two pieces "
                f"that share {run['overlap_pct']}% of the query")
        truth.update({
            "truth": description + _background_note(run),
            "true_copies_per_genome": copies,
            "true_left_contexts": copies, "true_right_contexts": copies,
            "variant_sites_1based": (f"copy2:{len(changes)} substitutions" if changes else ""),
        })
        return {
            "components": [("strain", genome, 1.0)], "bait": reference,
            "alleles": [(f"copy{i + 1}", s) for i, s in enumerate(copy_sequences)],
            "truth": truth,
        }

    # ---- test 4 ----------------------------------------------------------- #
    copies, layout, divergent = run["true_copies"], run["layout"], run["divergent"]
    reference = build_allele(length, rng)
    copy_sequences = [reference]
    divergence_notes = []
    for index in range(1, copies):
        if divergent:
            positions = rng.choice(length, size=round(T4_DIVERGENCE * length), replace=False)
            other, changes = substitute(reference, positions, rng)
            divergence_notes.append(f"copy{index + 1}:{len(changes)} substitutions")
        else:
            other = reference
        copy_sequences.append(other)
    alleles = [(f"copy{i + 1}", s) for i, s in enumerate(copy_sequences)]
    left_contexts = right_contexts = copies
    second_genome = None
    note = ""

    if layout == "dispersed":
        parts = _segments(copies + 1, rng, run)
        genome = parts[0]
        for index in range(copies):
            genome += copy_sequences[index] + parts[index + 1]
    elif layout.startswith("sharedflank"):
        shared_bp = int(layout.replace("sharedflank", ""))
        parts = _segments(3, rng, run)
        shared = parts[0][-shared_bp:]
        genome = (parts[0] + copy_sequences[0] + parts[1] + shared
                  + copy_sequences[1] + parts[2])
        note = (f"the {shared_bp} bp upstream of both copies are identical; "
                "the flanks differ beyond them and on the other side")
    elif layout == "tandem":
        parts = _segments(2, rng, run)
        genome = parts[0] + copy_sequences[0] + copy_sequences[1] + parts[1]
        left_contexts = right_contexts = 1
        note = "direct tandem duplicate without spacer; one external context per side"
    elif layout == "flankrepeat":
        parts = _segments(3, rng, run)
        repeat = random_dna(T4_FLANK_REPEAT_BP, rng)
        spacer = random_dna(T4_FLANK_REPEAT_SPACER, rng)
        genome = (parts[0] + repeat + parts[1] + repeat + spacer
                  + copy_sequences[0] + parts[2])
        note = (f"a {T4_FLANK_REPEAT_BP}-bp repeat lies {T4_FLANK_REPEAT_SPACER} bp "
                "upstream of the locus and once more elsewhere in the genome")
    elif layout == "strainmix":
        parts = _segments(copies + 1, rng, run)
        genome = parts[0]
        for index in range(copies):
            genome += copy_sequences[index] + parts[index + 1]
        site = len(parts[0]) - T4_FLANK_SITE_DISTANCE
        original = genome[site]
        alternative = [b for b in DNA if b != original][int(rng.integers(0, 3))]
        second_genome = genome[:site] + alternative + genome[site + 1:]
        note = (f"50:50 mixture of two strains that differ at one site "
                f"{T4_FLANK_SITE_DISTANCE} bp upstream of copy 1 "
                f"(genome position {site + 1}, {original}>{alternative}); "
                "each strain has the stated number of copies")
    else:
        raise ValueError(layout)

    if second_genome is None:
        components = [("strain", genome, 1.0)]
    else:
        components = [("strain_1", genome, 0.5), ("strain_2", second_genome, 0.5)]
    truth.update({
        "truth": (f"{copies} cop{'y' if copies == 1 else 'ies'}, layout {layout}, "
                  f"copies {'1% divergent' if divergent else 'identical'}"
                  + (f"; {note}" if note else "") + _background_note(run)),
        "true_copies_per_genome": copies,
        "true_left_contexts": left_contexts, "true_right_contexts": right_contexts,
        "variant_sites_1based": ";".join(divergence_notes),
    })
    return {"components": components, "bait": reference, "alleles": alleles, "truth": truth}


def _test3_allele(scenario: str, reference: str, rng):
    """Sample allele of test 3, its description and the events applied."""
    length = len(reference)
    margin = T3_INTERIOR_MARGIN
    if scenario == "control_same_as_bait":
        return reference, "allele identical to the bait", ""
    if scenario.startswith("sub_ends_"):
        distance = int(scenario.split("_")[2].replace("bp", ""))
        five = int(rng.integers(0, distance))
        three = length - 1 - int(rng.integers(0, distance))
        sample, changes = substitute(reference, [five, three], rng)
        return (
            sample,
            f"one substitution within {distance} bp of each end of the bait",
            ";".join(f"{p}{r}>{a}" for p, r, a in changes),
        )
    first_codon, last_codon = margin // 3 + 1, (length - margin) // 3 - 1
    if scenario == "del_3bp":
        codon = int(rng.integers(first_codon, last_codon))
        sample = reference[:codon * 3] + reference[codon * 3 + 3:]
        return sample, "in-frame deletion of one codon", f"del{codon * 3 + 1}-{codon * 3 + 3}"
    if scenario == "ins_12bp":
        codon = int(rng.integers(first_codon, last_codon))
        inserted = "".join(SENSE_CODONS[i] for i in rng.integers(0, len(SENSE_CODONS), 4))
        sample = reference[:codon * 3] + inserted + reference[codon * 3:]
        return sample, "in-frame insertion of four codons", f"ins{codon * 3}^{inserted}"
    if scenario == "del_1bp":
        position = int(rng.integers(margin, length - margin))
        sample = reference[:position] + reference[position + 1:]
        return sample, "deletion of one base (frameshift)", f"del{position + 1}"
    if scenario == "compensating_indels_30bp":
        position = int(rng.integers(margin, length - margin - 30))
        inserted = DNA[int(rng.integers(0, 4))]
        sample = (reference[:position] + reference[position + 1:position + 30]
                  + inserted + reference[position + 30:])
        return (
            sample,
            "deletion of one base and insertion of one base 30 bp downstream "
            "(length kept, ten codons altered)",
            f"del{position + 1};ins{position + 30}^{inserted}",
        )
    if scenario == "premature_stop":
        candidates = []
        for codon in range(first_codon, last_codon):
            triplet = reference[codon * 3:codon * 3 + 3]
            for offset in range(3):
                for base in DNA:
                    if base == triplet[offset]:
                        continue
                    if triplet[:offset] + base + triplet[offset + 1:] in STOP_CODONS:
                        candidates.append((codon * 3 + offset, base))
        position, base = candidates[int(rng.integers(0, len(candidates)))]
        sample = reference[:position] + base + reference[position + 1:]
        return (
            sample, "one substitution that creates a stop codon in frame +1",
            f"{position + 1}{reference[position]}>{base}",
        )
    raise ValueError(scenario)


# --------------------------------------------------------------------------- #
# External programs
# --------------------------------------------------------------------------- #
def _run(command, log_path: Path, cwd: Path | None = None) -> int:
    with log_path.open("a") as log:
        log.write("$ " + " ".join(str(c) for c in command) + "\n")
        log.flush()
        completed = subprocess.run(
            [str(c) for c in command], cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
        )
        log.write(f"[exit {completed.returncode}]\n")
    return completed.returncode


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, chunks = [], None, []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records


def assembly_stats(records) -> dict:
    lengths = sorted((len(s) for _n, s in records), reverse=True)
    total = sum(lengths)
    running, n50 = 0, 0
    for value in lengths:
        running += value
        if running * 2 >= total:
            n50 = value
            break
    return {
        "draft_contigs": len(lengths), "draft_total_bp": total, "draft_n50": n50,
        "draft_longest_bp": lengths[0] if lengths else 0,
    }


def assemble_draft(r1: Path, r2: Path, out_dir: Path, threads: int, memory_gb: int,
                   log: Path) -> tuple[int, Path, Path]:
    """SPAdes ``--isolate`` draft, as ``workflows/run_low_copy_23S.py``."""
    command = [
        "spades.py", "--isolate", "-1", r1, "-2", r2, "-o", out_dir,
        "--threads", threads, "--memory", memory_gb,
    ]
    code = _run(command, log)
    return (code, out_dir / "scaffolds.fasta",
            out_dir / "assembly_graph_after_simplification.gfa")


def prune_spades(directory: Path) -> None:
    """Drop the k-mer working directories; keep contigs, graphs and logs."""
    if not directory.is_dir():
        return
    for child in directory.iterdir():
        if child.is_dir() and (re.fullmatch(r"K\d+", child.name)
                               or child.name in {"tmp", "corrected", "misc", "pipeline_state"}):
            shutil.rmtree(child, ignore_errors=True)


def discover(draft: Path, bait: Path, bait_length: int, work: Path, log: Path):
    """Discovery intervals for the copy-number commands.

    The rule of ``workflows/run_low_copy_23S.py::_discover``: BLASTN of the bait
    against the draft at >= 80% identity, hits of at least max(200 bp, 10% of
    the bait), the best-scoring hit kept wherever hits overlap on the subject.
    Returns regions (``contig:start-end``), query spans and the hit rows.
    """
    database = work / "draft_db"
    _run(["makeblastdb", "-in", draft, "-dbtype", "nucl", "-out", database], log)
    output = work / "bait_vs_draft.tsv"
    fields = "qseqid sseqid pident length qstart qend sstart send bitscore"
    _run([
        "blastn", "-query", bait, "-db", database, "-out", output,
        "-outfmt", f"6 {fields}", "-max_target_seqs", "20", "-perc_identity", "80",
    ], log)
    hits = []
    if output.exists():
        for line in output.read_text().splitlines():
            _q, sid, pid, hit_len, qs, qe, ss, se, bitscore = line.split("\t")
            if int(hit_len) < max(200, int(0.10 * bait_length)):
                continue
            hits.append({
                "contig": sid, "subject_start": min(int(ss), int(se)),
                "subject_end": max(int(ss), int(se)),
                "strand": "+" if int(ss) <= int(se) else "-",
                "query_start": min(int(qs), int(qe)), "query_end": max(int(qs), int(qe)),
                "bitscore": float(bitscore), "identity": float(pid),
            })
    retained: list[dict] = []
    for hit in sorted(hits, key=lambda item: item["bitscore"], reverse=True):
        if not any(old["contig"] == hit["contig"]
                   and hit["subject_start"] <= old["subject_end"]
                   and old["subject_start"] <= hit["subject_end"] for old in retained):
            retained.append(hit)
    retained.sort(key=lambda item: (item["contig"], item["subject_start"]))
    regions = [f"{h['contig']}:{h['subject_start']}-{h['subject_end']}" for h in retained]
    spans = [(h["query_start"], h["query_end"]) for h in retained]
    return regions, spans, retained


def cut_draft(draft: Path, bait: Path, bait_length: int, overlap_pct: int, work: Path,
              log: Path) -> tuple[Path | None, dict]:
    """Cut the contig that holds the whole locus into two overlapping pieces.

    Piece A ends at query position ``b`` and piece B starts at query position
    ``a``; they share ``b - a + 1 = round(overlap * L)`` query bases, centred
    on the query.  At 0% the pieces are adjacent.  Every other contig is kept.
    """
    cut_work = work / "uncut_discovery"
    cut_work.mkdir(exist_ok=True)
    _regions, _spans, hits = discover(draft, bait, bait_length, cut_work, log)
    full = [h for h in hits if h["query_start"] == 1 and h["query_end"] == bait_length]
    info = {"uncut_discovery_hits": len(hits), "uncut_full_length_hits": len(full)}
    if len(full) != 1 or len(hits) != 1:
        info["cut_status"] = "NOT_CUT_NO_SINGLE_FULL_LENGTH_HIT"
        return None, info
    hit = full[0]
    overlap_bp = round(overlap_pct / 100.0 * bait_length)
    a = (bait_length - overlap_bp) // 2 + 1
    b = a + overlap_bp - 1
    records = read_fasta(draft)
    out_records = []
    for name, sequence in records:
        if name != hit["contig"]:
            out_records.append((name, sequence))
            continue
        if hit["strand"] == "-":
            sequence = reverse_complement(sequence)
            start0 = len(sequence) - hit["subject_end"]
        else:
            start0 = hit["subject_start"] - 1
        piece_a = sequence[:start0 + b]
        piece_b = sequence[start0 + a - 1:]
        out_records.append((f"{name}_cutA", piece_a))
        out_records.append((f"{name}_cutB", piece_b))
        info.update({
            "cut_status": "CUT", "cut_contig": name, "cut_contig_bp": len(sequence),
            "cut_query_a": a, "cut_query_b": b, "cut_overlap_bp": overlap_bp,
            "cut_piece_a_bp": len(piece_a), "cut_piece_b_bp": len(piece_b),
        })
    cut_path = work / "draft_cut.fasta"
    write_fasta(cut_path, out_records)
    return cut_path, info


def cut_draft_divergent(draft: Path, bait: Path, bait_length: int, overlap_pct: int,
                        work: Path, log: Path) -> tuple[Path | None, dict]:
    """Keep copy 1 over query 1..b and copy 2 over query a..L; remove the rest.

    Needs a draft in which both copies are complete: the copy identical to the
    bait (copy 1) and the copy that differs from it (copy 2).  The bases of
    copy 1 after ``b`` and of copy 2 before ``a`` are excised and the contig is
    split there.  ``a`` and ``b`` are those of ``cut_draft``.
    """
    cut_work = work / "uncut_discovery"
    cut_work.mkdir(exist_ok=True)
    _regions, _spans, hits = discover(draft, bait, bait_length, cut_work, log)
    # A substitution close to an end of copy 2 shortens the local alignment by a
    # few bases.  The copies differ by substitutions only, so such a hit is
    # extended to the whole query, base for base.
    tolerance = 30
    full = []
    for hit in hits:
        missing_start = hit["query_start"] - 1
        missing_end = bait_length - hit["query_end"]
        if missing_start > tolerance or missing_end > tolerance:
            continue
        if hit["subject_end"] - hit["subject_start"] != hit["query_end"] - hit["query_start"]:
            continue
        extended = dict(hit)
        if hit["strand"] == "+":
            extended["subject_start"] -= missing_start
            extended["subject_end"] += missing_end
        else:
            extended["subject_start"] -= missing_end
            extended["subject_end"] += missing_start
        extended["query_start"], extended["query_end"] = 1, bait_length
        full.append(extended)
    exact = [h for h in full if h["identity"] == 100.0]
    info = {"uncut_discovery_hits": len(hits), "uncut_full_length_hits": len(full)}
    if len(hits) != 2 or len(full) != 2 or len(exact) != 1:
        info["cut_status"] = "NOT_CUT_NEEDS_TWO_COMPLETE_COPIES"
        return None, info
    copy_1 = exact[0]
    copy_2 = next(h for h in full if h is not copy_1)
    overlap_bp = round(overlap_pct / 100.0 * bait_length)
    a = (bait_length - overlap_bp) // 2 + 1
    b = a + overlap_bp - 1

    def subject_interval(hit, query_from, query_to):
        if query_from > query_to:
            return None
        if hit["strand"] == "+":
            return (hit["subject_start"] + query_from - 1, hit["subject_start"] + query_to - 1)
        return (hit["subject_end"] - (query_to - 1), hit["subject_end"] - (query_from - 1))

    removals: dict[str, list[tuple[int, int]]] = {}
    for hit, interval in ((copy_1, subject_interval(copy_1, b + 1, bait_length)),
                          (copy_2, subject_interval(copy_2, 1, a - 1))):
        if interval is not None:
            removals.setdefault(hit["contig"], []).append(interval)
    out_records = []
    pieces = 0
    for name, sequence in read_fasta(draft):
        if name not in removals:
            out_records.append((name, sequence))
            continue
        position = 0
        part = 0
        for start, end in sorted(removals[name]):
            if start - 1 > position:
                part += 1
                out_records.append((f"{name}_part{part}", sequence[position:start - 1]))
            position = end
        if position < len(sequence):
            part += 1
            out_records.append((f"{name}_part{part}", sequence[position:]))
        pieces += part
    info.update({
        "cut_status": "CUT", "cut_contig": copy_1["contig"], "cut_query_a": a,
        "cut_query_b": b, "cut_overlap_bp": overlap_bp, "cut_pieces": pieces,
        "cut_copy2_identity": copy_2["identity"],
    })
    cut_path = work / "draft_cut.fasta"
    write_fasta(cut_path, out_records)
    return cut_path, info


def run_locus_recon(sample: str, draft: Path, r1: Path, r2: Path, bait: Path, out_dir: Path,
                    threads: int, memory_gb: int, log: Path) -> tuple[int, dict]:
    """The standard workflow with default options (threads and memory only)."""
    sheet = out_dir.parent / "samples.tsv"
    sheet.write_text(
        "sample_id\tassembly_path\tr1_path\tr2_path\n"
        f"{sample}\t{draft}\t{r1}\t{r2}\n"
    )
    command = [
        "locus-recon", "--samplesheet", sheet, "--main-output-dir", out_dir,
        "--bait", bait, "--locus", LOCUS, "--threads", threads,
        "--memory-per-sample", memory_gb,
    ]
    code = _run(command, log)
    report = out_dir / f"locus_recon_report_{LOCUS}.tsv"
    row: dict = {}
    if report.exists():
        with report.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if rows:
            row = rows[0]
    return code, row


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #
def compare_sequences(returned: str, truth: str) -> dict:
    """Identity of a returned sequence to one truth sequence."""
    result = {
        "identical": False, "orientation": "+", "incomplete_identical": False,
        "mismatches": None, "inserted_bp": None, "deleted_bp": None,
        "missing_5prime_bp": None, "missing_3prime_bp": None,
    }
    if not returned:
        return result
    candidate = returned
    if returned != truth and reverse_complement(returned) == truth:
        candidate, result["orientation"] = reverse_complement(returned), "-"
    if candidate == truth:
        result.update({"identical": True, "mismatches": 0, "inserted_bp": 0,
                       "deleted_bp": 0, "missing_5prime_bp": 0, "missing_3prime_bp": 0})
        return result
    position = truth.find(candidate)
    if position >= 0:
        result.update({
            "incomplete_identical": True, "mismatches": 0, "inserted_bp": 0,
            "deleted_bp": len(truth) - len(candidate), "missing_5prime_bp": position,
            "missing_3prime_bp": len(truth) - len(candidate) - position,
        })
        return result
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 1, -1
    aligner.open_gap_score, aligner.extend_gap_score = -2, -1
    alignment = aligner.align(truth, candidate)[0]
    mismatches = 0
    truth_aligned = candidate_aligned = 0
    for (t0, t1), (c0, c1) in zip(*alignment.aligned):
        truth_aligned += t1 - t0
        candidate_aligned += c1 - c0
        mismatches += sum(1 for x, y in zip(truth[t0:t1], candidate[c0:c1]) if x != y)
    result.update({
        "mismatches": int(mismatches),
        "deleted_bp": int(len(truth) - truth_aligned),       # truth bases absent
        "inserted_bp": int(len(candidate) - candidate_aligned),  # extra bases
    })
    return result


def _plain(value):
    """Numpy scalars and containers as plain JSON types."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def flag_names(flags: str) -> list[str]:
    return FLAG_NAME.findall(flags or "")


REPORT_FIELDS = (
    "status", "workflow_status", "result_disposition", "qc_confidence",
    "allele_length", "nearest_allele", "exact_known_allele", "best_identity_pct",
    "validation_coverage_pct", "contigs_hit", "target_regions", "spades_mode",
    "qc_length_delta", "qc_internal_stops", "qc_coding_frame_used", "qc_length_mod3",
    "remap_mean_depth", "remap_breadth_pct", "uncertain_base_count",
    "uncertain_base_fraction", "interior_uncertain_base_count",
    "interior_uncertain_base_fraction", "internal_zero_depth_positions",
    "terminal_margin_bp", "span_clipped_bp", "strand_balance_pct",
    "candidate_mixed_sites", "mixed_site_count", "strand_biased_sites",
    "reference_discordant_sites", "max_alt_fraction", "median_mixed_fraction",
    "mixture_detected", "bitscore_margin", "qc_flags", "message",
)

DEPTH_FIELDS = (
    "locus_regions", "region_gene_spans", "region_ratios", "assembly_region_count",
    "locus_median_depth", "locus_median_depth_unfiltered", "backbone_median_depth",
    "ratio", "ratio_ci_low", "ratio_ci_high", "depth_call", "ambiguity_index", "reliable",
    "dosage_estimate", "dosage_ci_low", "dosage_ci_high", "dosage_lower", "dosage_upper",
    "dosage_status", "dosage_flags", "geometry_expected_dosage", "gene_overlap_bp",
    "gene_overlap_fraction", "query_coverage_fraction", "query_multiplicity_max",
    "bootstrap_block_size", "bootstrap_locus_blocks", "bootstrap_backbone_blocks_used",
)
GRAPH_FIELDS = (
    "graph_target_segments", "graph_bait_coverage_fraction", "left_context_count",
    "right_context_count", "graph_context_count", "graph_context_lower_bound",
    "graph_context_status", "graph_context_flags", "graph_min_context_bp",
)
CONSENSUS_FIELDS = (
    "copy_number_call", "copy_number_kind", "copy_number_method", "copy_number_status",
    "copy_number_lower_bound", "copy_number_upper_bound", "consensus_flags",
)


def _flat(value):
    if isinstance(value, (list, tuple)):
        if value and isinstance(value[0], (list, tuple)):
            return ";".join(f"{a}-{b}" for a, b in value)
        return ";".join(str(item) for item in value)
    return value


def copy_number_block(sample_dir: Path, lr_dir: Path, sample: str, draft: Path, gfa: Path,
                      bait: Path, bait_length: int, true_copies: int, with_no_graph: bool,
                      log: Path) -> tuple[dict, list[dict]]:
    """Depth ratio, copy-number consensus and the block-size series of one run."""
    out: dict = {}
    block_rows: list[dict] = []
    work = sample_dir / "copy_number"
    work.mkdir(exist_ok=True)
    regions, spans, hits = discover(draft, bait, bait_length, work, log)
    out["cn_discovery_intervals"] = len(regions)
    out["cn_discovery_identities"] = ";".join(f"{h['identity']:.2f}" for h in hits)
    if not regions:
        out["cn_status"] = "NO_DISCOVERY_INTERVAL"
        return out, block_rows
    bam = lr_dir / sample / "full_map.sorted.bam"
    if not bam.exists():
        out["cn_status"] = "NO_FULL_MAP_BAM"
        return out, block_rows
    depth_args: list = ["--bam", bam]
    for region, (start, end) in zip(regions, spans):
        depth_args += ["--locus", region, "--gene-span", f"{start}-{end}"]
    depth_args += ["--gene-length", bait_length, "--reference", draft]

    depth_tsv = work / "depth_ratio.tsv"
    code = _run(["locus-recon-depth-ratio", *depth_args, "--tsv", depth_tsv], log)
    out["depth_ratio_exit"] = code
    if depth_tsv.exists():
        with depth_tsv.open(newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        for key in DEPTH_FIELDS:
            out[f"dr_{key}"] = row.get(key, "")

    def consensus(label: str, extra: list) -> None:
        json_path = work / f"copy_number_{label}.json"
        code = _run([
            "locus-recon-copy-number", *depth_args, *extra,
            "--output-json", json_path, "--output-tsv", work / f"copy_number_{label}.tsv",
        ], log)
        out[f"{label}_exit"] = code
        if not json_path.exists():
            return
        record = json.loads(json_path.read_text())
        if label == "cn":
            for key in DEPTH_FIELDS:
                out[f"cn_{key}"] = _flat(record.get(key))
            for key in GRAPH_FIELDS:
                out[key] = _flat(record.get(key))
        for key in CONSENSUS_FIELDS:
            out[f"{label}_{key}"] = _flat(record.get(key))

    consensus("cn", ["--gfa", gfa, "--bait", bait])
    if with_no_graph:
        consensus("cn_nograph", [])
    out["cn_status"] = "DONE"

    from locus_recon.depth_ratio import estimate_locus_copy_number
    for block in BLOCK_SIZES:
        result = estimate_locus_copy_number(
            bam, regions, gene_spans=spans, gene_length=bait_length,
            block=block, n_boot=BOOTSTRAP_RESAMPLES, seed=BOOTSTRAP_SEED,
        )
        block_rows.append({
            "block_bp": block,
            "true_copies": true_copies,
            "assembly_region_count": result.assembly_region_count,
            "ratio": result.ratio, "ratio_ci_low": result.ratio_ci_low,
            "ratio_ci_high": result.ratio_ci_high,
            "ratio_ci_width": round(result.ratio_ci_high - result.ratio_ci_low, 3),
            "depth_call": result.depth_call,
            "dosage_estimate": result.dosage_estimate,
            "dosage_ci_low": result.dosage_ci_low, "dosage_ci_high": result.dosage_ci_high,
            "dosage_ci_width": (
                round(result.dosage_ci_high - result.dosage_ci_low, 3)
                if result.dosage_ci_low is not None else ""),
            "dosage_status": result.dosage_status,
            "dosage_ci_contains_true_copies": (
                int(result.dosage_ci_low <= true_copies <= result.dosage_ci_high)
                if result.dosage_ci_low is not None else ""),
            "ratio_ci_contains_true_copies": int(
                result.ratio_ci_low <= true_copies <= result.ratio_ci_high),
            "bootstrap_locus_blocks": result.bootstrap_locus_blocks,
            "bootstrap_backbone_blocks_total": result.bootstrap_backbone_blocks_total,
            "bootstrap_backbone_blocks_used": result.bootstrap_backbone_blocks_used,
        })
    return out, block_rows


# --------------------------------------------------------------------------- #
# One run
# --------------------------------------------------------------------------- #
def process_run(run: dict, work_dir: str, threads: int, memory_gb: int) -> dict:
    run_id = run["run_id"]
    sample_dir = Path(work_dir) / "runs" / f"test{run['test']}" / run_id
    done = sample_dir / "result.json"
    if done.exists():
        return json.loads(done.read_text())
    if sample_dir.exists():
        shutil.rmtree(sample_dir)
    sample_dir.mkdir(parents=True)
    log = sample_dir / "commands.log"
    started = time.time()
    record: dict = OrderedDict(run)
    record["seed_genome"] = seed_for(run_id, "genome")
    record["seed_reads"] = seed_for(run_id, "reads")
    checksums: list[dict] = []
    try:
        case = build_case(run)
        record.update(case["truth"])
        truth_dir = sample_dir / "truth"
        truth_dir.mkdir()
        bait = truth_dir / "bait.fasta"
        write_fasta(bait, [(BAIT_ID, case["bait"])])
        write_fasta(truth_dir / "truth_alleles.fasta",
                    [(f"{run_id}|{name}", seq) for name, seq in case["alleles"]])
        write_fasta(truth_dir / "truth_genomes.fasta",
                    [(f"{run_id}|{name}|fraction={fraction}", genome)
                     for name, genome, fraction in case["components"]])
        record["genome_bp"] = len(case["components"][0][1])

        r1, r2 = sample_dir / "reads_R1.fastq.gz", sample_dir / "reads_R2.fastq.gz"
        rng = np.random.default_rng(record["seed_reads"])
        read_1, read_2 = simulate_sample(
            [(genome, fraction) for _n, genome, fraction in case["components"]],
            run["depth"], rng)
        write_fastq_pair(r1, r2, read_1, read_2, run_id)
        record["read_pairs"] = int(len(read_1))
        del read_1, read_2

        code, draft, gfa = assemble_draft(
            r1, r2, sample_dir / "draft_spades", threads, memory_gb, log)
        record["draft_spades_exit"] = code
        prune_spades(sample_dir / "draft_spades")
        if code != 0 or not draft.exists() or draft.stat().st_size == 0:
            record["run_status"] = "DRAFT_ASSEMBLY_FAILED"
            raise _Stop()
        record.update(assembly_stats(read_fasta(draft)))
        draft_used = draft
        if run["test"] == 5:
            cut_source = draft
            record["draft_source"] = "SPAdes scaffolds, one contig cut in two"
            cutter = cut_draft
            if run.get("divergent"):
                # SPAdes does not always return both divergent copies complete,
                # so this added arm starts from the true chromosome, as the
                # constructed drafts of generate_mock_dataset.py do.  The SPAdes
                # graph of the same reads is still the graph given to
                # locus-recon-copy-number.
                cut_source = sample_dir / "truth_draft.fasta"
                write_fasta(cut_source, [("truth_chromosome", case["components"][0][1])])
                record["draft_source"] = "true chromosome with the stated parts removed"
                cutter = cut_draft_divergent
            cut, info = cutter(cut_source, bait, len(case["bait"]), run["overlap_pct"],
                               sample_dir, log)
            record.update(info)
            if cut is None:
                record["run_status"] = info["cut_status"]
                raise _Stop()
            draft_used = cut

        lr_dir = sample_dir / "locus_recon"
        code, row = run_locus_recon(run_id, draft_used, r1, r2, bait, lr_dir,
                                    threads, memory_gb, log)
        record["locus_recon_exit"] = code
        for key in REPORT_FIELDS:
            record[f"lr_{key}"] = row.get(key, "")
        names = flag_names(row.get("qc_flags", ""))
        record["lr_flag_names"] = ";".join(names)
        record["tool_version"] = row.get("tool_version", "")

        returned = ""
        allele_file = row.get("allele_file", "")
        if allele_file and (lr_dir / allele_file).exists():
            fasta = read_fasta(lr_dir / allele_file)
            if fasta:
                returned = fasta[0][1]
        record["sequence_returned"] = int(bool(returned))
        record["returned_length"] = len(returned)
        if returned:
            write_fasta(sample_dir / "returned.fasta", [(run_id, returned)])
        best = None
        for name, truth_sequence in case["alleles"]:
            comparison = compare_sequences(returned, truth_sequence)
            comparison["allele"] = name
            key = (comparison["identical"], comparison["incomplete_identical"],
                   -((comparison["mismatches"] or 0) + (comparison["inserted_bp"] or 0)
                     + (comparison["deleted_bp"] or 0)) if returned else -10**9)
            if best is None or key > best[0]:
                best = (key, comparison)
        comparison = best[1]
        record["cmp_each_allele"] = ";".join(
            "{}:identical={},mismatches={},inserted={},deleted={}".format(
                name, int(c["identical"]), c["mismatches"], c["inserted_bp"], c["deleted_bp"])
            for name, c in ((n, compare_sequences(returned, s)) for n, s in case["alleles"])
        ) if returned else ""
        record["exact"] = int(comparison["identical"])
        record["matches_allele"] = comparison["allele"] if comparison["identical"] else ""
        record["closest_allele"] = comparison["allele"] if returned else ""
        for key in ("orientation", "incomplete_identical", "mismatches", "inserted_bp",
                    "deleted_bp", "missing_5prime_bp", "missing_3prime_bp"):
            value = comparison[key]
            record[f"cmp_{key}"] = int(value) if isinstance(value, bool) else value

        if run["test"] in (4, 5):
            extra, block_rows = copy_number_block(
                sample_dir, lr_dir, run_id, draft_used, gfa, bait, len(case["bait"]),
                run["true_copies"], with_no_graph=True, log=log)
            record.update(extra)
            record["block_rows"] = block_rows
        record["run_status"] = "DONE"

        bam = lr_dir / run_id / "full_map.sorted.bam"
        for role, path in (("reads_R1", r1), ("reads_R2", r2), ("bait", bait),
                           ("truth_alleles", truth_dir / "truth_alleles.fasta"),
                           ("truth_genomes", truth_dir / "truth_genomes.fasta"),
                           ("draft_scaffolds", draft), ("draft_gfa", gfa),
                           ("draft_used", draft_used), ("full_map_bam", bam)):
            if Path(path).exists():
                checksums.append({"run_id": run_id, "role": role, "file": Path(path).name,
                                  "bytes": Path(path).stat().st_size,
                                  "sha256": sha256_file(path)})
        prune_spades(lr_dir / run_id / "spades_local")
    except _Stop:
        pass
    except Exception as exc:  # recorded, never hidden
        record["run_status"] = "RUNNER_ERROR"
        record["runner_error"] = f"{type(exc).__name__}: {exc}"
        (sample_dir / "runner_error.txt").write_text(traceback.format_exc())
    record["checksums"] = checksums
    record["elapsed_sec"] = round(time.time() - started, 1)
    record["remote_dir"] = str(sample_dir)
    record = _plain(record)
    write_json(done, record)
    return record


class _Stop(Exception):
    """A run that ended at a recorded state before the workflow finished."""


# --------------------------------------------------------------------------- #
# Tables
# --------------------------------------------------------------------------- #
COMMON = [
    "run_id", "test", "scenario", "replicate", "depth", "locus_bp", "truth",
    "variant_sites_1based", "true_copies_per_genome", "seed_genome", "seed_reads",
    "genome_bp", "read_pairs", "run_status", "draft_spades_exit", "draft_contigs",
    "draft_total_bp", "draft_n50", "locus_recon_exit",
]
LR = ["lr_" + f for f in REPORT_FIELDS if f not in ("qc_flags", "message")] + [
    "lr_flag_names", "lr_qc_flags", "sequence_returned", "returned_length", "exact",
    "matches_allele", "closest_allele", "cmp_orientation", "cmp_incomplete_identical",
    "cmp_mismatches", "cmp_inserted_bp", "cmp_deleted_bp", "cmp_missing_5prime_bp",
    "cmp_missing_3prime_bp", "cmp_each_allele",
]
CN = (["cn_status", "cn_discovery_intervals", "cn_discovery_identities"]
      + ["dr_" + f for f in DEPTH_FIELDS] + list(GRAPH_FIELDS)
      + ["cn_" + f for f in CONSENSUS_FIELDS]
      + ["cn_nograph_" + f for f in CONSENSUS_FIELDS]
      + ["depth_ratio_exit", "cn_exit", "cn_nograph_exit"])
TABLES = {
    1: ("test1_mixtures_few_sites.tsv",
        COMMON[:5] + ["sites", "minor_fraction_pct"] + COMMON[5:] + LR),
    2: ("test2_strand_bias_pure.tsv", COMMON[:5] + ["nominal_locus_bp"] + COMMON[5:] + LR),
    3: ("test3_ends_indels_pseudogenes.tsv",
        COMMON + ["truth_length", "truth_length_mod3", "bait_length_mod3",
                  "truth_internal_stops_frame_plus1", "truth_internal_stops_min_six_frames",
                  "truth_frame_intact"] + LR),
    4: ("test4_graph_contexts.tsv",
        COMMON[:5] + ["background", "true_copies", "layout", "divergent", "true_left_contexts",
                      "true_right_contexts"] + COMMON[5:] + CN + LR),
    5: ("test5_overlapping_spans.tsv",
        COMMON[:5] + ["arm", "background", "true_copies", "divergent", "overlap_pct",
                      "draft_source", "cut_status", "cut_contig_bp",
                      "cut_query_a", "cut_query_b", "cut_overlap_bp", "cut_piece_a_bp",
                      "cut_piece_b_bp", "cut_pieces", "cut_copy2_identity",
                      "uncut_discovery_hits"] + COMMON[5:] + CN + LR),
}


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: ("" if row.get(c) is None else row.get(c)) for c in columns})


def _has(row: dict, flag: str) -> bool:
    return flag in (row.get("lr_flag_names") or "").split(";")


def _num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _range(values) -> str:
    values = [v for v in values if v is not None]
    if not values:
        return ""
    return f"{statistics.median(values):.3f} ({min(values):.3f}-{max(values):.3f})"


def _counts(values) -> str:
    counter = Counter("" if v is None else str(v) for v in values)
    return "; ".join(f"{key or 'none'}={count}" for key, count in sorted(counter.items()))


DESCRIPTIVE = (
    "PATCHY_REMAP_SUPPORT", "MIXED_SITE_BELOW_MIXTURE_CALL", "FRAME_LENGTH_SHIFT_TRUNCATED",
    "PLACEHOLDER_FRAMESHIFT_EXPLAINS_STOPS", "PLACEHOLDER_JUNCTION_",
    "MODEST_REMAP_PATCHINESS", "REMAP_PATCHINESS_BELOW_HIGH", "CATALOGUE_",
    "POSSIBLE_NOVEL_ALLELE",
)


def aggregate(records: list[dict], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: list[dict] = []

    def add(test, cell, metric, value, denominator, source):
        summary.append({"test": test, "cell": cell, "metric": metric, "value": value,
                        "denominator": denominator, "source_file": source})

    for record in records:
        record["lr_qc_flags"] = record.get("lr_qc_flags", "")
    by_test = {t: [r for r in records if r["test"] == t] for t in (1, 2, 3, 4, 5)}

    # manifest: every run with its truth
    manifest_columns = [
        "run_id", "test", "scenario", "replicate", "depth", "locus_bp", "sites",
        "minor_fraction_pct", "nominal_locus_bp", "background", "arm", "true_copies", "layout",
        "divergent",
        "overlap_pct", "truth", "variant_sites_1based", "true_copies_per_genome",
        "true_left_contexts", "true_right_contexts", "seed_genome", "seed_reads",
        "genome_bp", "read_pairs", "run_status", "elapsed_sec", "remote_dir",
    ]
    _write_tsv(out_dir / "targeted_manifest.tsv", records, manifest_columns)
    for test, (name, columns) in TABLES.items():
        _write_tsv(out_dir / name, by_test[test], columns)

    checksum_rows = [c for r in records for c in r.get("checksums", [])]
    _write_tsv(out_dir / "input_checksums.tsv", checksum_rows,
               ["run_id", "role", "file", "bytes", "sha256"])

    # ---- test 1 ----------------------------------------------------------- #
    source = TABLES[1][0]
    cells = []
    groups = OrderedDict()
    for r in by_test[1]:
        groups.setdefault((r["scenario"], r.get("sites"), r.get("minor_fraction_pct")), []).append(r)
    for (scenario, sites, fraction), rows in groups.items():
        accepted = [r for r in rows if r.get("lr_result_disposition") == "PASS"]
        held = [r for r in rows if r.get("lr_result_disposition") == "HOLD"]
        cell = {
            "scenario": scenario, "sites": sites, "minor_fraction_pct": fraction,
            "runs": len(rows),
            "sequence_returned": sum(int(r.get("sequence_returned") or 0) for r in rows),
            "held": len(held),
            "held_as_mixture": sum(1 for r in held if r.get("lr_mixture_detected") == "true"),
            "held_with_STRAND_BIASED_ALT_SUPPORT": sum(
                1 for r in held if _has(r, "STRAND_BIASED_ALT_SUPPORT")),
            "sent_to_review": sum(1 for r in rows if r.get("lr_result_disposition") == "REVIEW"),
            "accepted": len(accepted),
            "accepted_with_MIXED_SITE_BELOW_MIXTURE_CALL": sum(
                1 for r in accepted if _has(r, "MIXED_SITE_BELOW_MIXTURE_CALL")),
            "accepted_without_any_mixed_site_flag": sum(
                1 for r in accepted if not _has(r, "MIXED_SITE_BELOW_MIXTURE_CALL")),
            "exact_major_allele": sum(1 for r in rows if r.get("matches_allele") in ("major", "allele")),
            "exact_minor_allele": sum(1 for r in rows if r.get("matches_allele") == "minor"),
            "not_exact": sum(1 for r in rows if int(r.get("sequence_returned") or 0)
                             and not int(r.get("exact") or 0)),
            "accepted_not_exact": sum(1 for r in accepted if not int(r.get("exact") or 0)),
            "mixed_site_count_values": _counts(r.get("lr_mixed_site_count") for r in rows),
            "dispositions": _counts(r.get("lr_result_disposition") for r in rows),
        }
        cells.append(cell)
        label = (f"sites={sites};minor={fraction}%" if scenario == "mixture" else "pure_control")
        for metric in ("held_as_mixture", "held", "accepted",
                       "accepted_with_MIXED_SITE_BELOW_MIXTURE_CALL", "accepted_not_exact"):
            add(1, label, metric, cell[metric], len(rows), "test1_cells.tsv")
    _write_tsv(out_dir / "test1_cells.tsv", cells, list(cells[0]) if cells else [])
    limits = []
    for sites in T1_SITES:
        per_fraction = {c["minor_fraction_pct"]: c for c in cells
                        if c["scenario"] == "mixture" and c["sites"] == sites}
        all_held, lowest = None, None
        for fraction in sorted(per_fraction, reverse=True):
            cell = per_fraction[fraction]
            if cell["held_as_mixture"] == cell["runs"]:
                all_held = fraction
            else:
                break
        for fraction in sorted(per_fraction):
            if per_fraction[fraction]["held_as_mixture"] > 0:
                lowest = fraction
                break
        limits.append({
            "sites": sites,
            "lowest_minor_fraction_pct_with_any_replicate_held_as_mixture": lowest if lowest is not None else "none",
            "lowest_minor_fraction_pct_from_which_all_replicates_held_as_mixture": all_held if all_held is not None else "none",
            "held_as_mixture_by_fraction": "; ".join(
                f"{f}%={per_fraction[f]['held_as_mixture']}/{per_fraction[f]['runs']}"
                for f in sorted(per_fraction)),
            "accepted_by_fraction": "; ".join(
                f"{f}%={per_fraction[f]['accepted']}/{per_fraction[f]['runs']}"
                for f in sorted(per_fraction)),
        })
        add(1, f"sites={sites}", "detection_limit_all_replicates_held_as_mixture_pct",
            limits[-1]["lowest_minor_fraction_pct_from_which_all_replicates_held_as_mixture"],
            "", "test1_detection_limit.tsv")
    _write_tsv(out_dir / "test1_detection_limit.tsv", limits, list(limits[0]) if limits else [])

    # ---- test 2 ----------------------------------------------------------- #
    cells = []
    groups = OrderedDict()
    for r in by_test[2]:
        groups.setdefault((r.get("nominal_locus_bp"), r["locus_bp"], r["depth"]), []).append(r)
    for (nominal, length, depth), rows in groups.items():
        other = Counter()
        for r in rows:
            for name in (r.get("lr_flag_names") or "").split(";"):
                if name and name != "STRAND_BIASED_ALT_SUPPORT":
                    other[name] += 1
        tier_flag_runs = sum(
            1 for r in rows
            if any(n and not n.startswith(DESCRIPTIVE)
                   for n in (r.get("lr_flag_names") or "").split(";")))
        cell = {
            "nominal_locus_bp": nominal, "locus_bp": length, "depth": depth, "runs": len(rows),
            "sequence_returned": sum(int(r.get("sequence_returned") or 0) for r in rows),
            "exact": sum(int(r.get("exact") or 0) for r in rows),
            "runs_with_STRAND_BIASED_ALT_SUPPORT": sum(
                1 for r in rows if _has(r, "STRAND_BIASED_ALT_SUPPORT")),
            "runs_with_any_tier_affecting_flag": tier_flag_runs,
            "runs_with_any_flag": sum(1 for r in rows if r.get("lr_flag_names")),
            "runs_with_candidate_mixed_site": sum(
                1 for r in rows if (_num(r.get("lr_candidate_mixed_sites")) or 0) > 0),
            "accepted": sum(1 for r in rows if r.get("lr_result_disposition") == "PASS"),
            "sent_to_review": sum(1 for r in rows if r.get("lr_result_disposition") == "REVIEW"),
            "held": sum(1 for r in rows if r.get("lr_result_disposition") == "HOLD"),
            "other_flags": "; ".join(f"{k}={v}" for k, v in sorted(other.items())),
        }
        cells.append(cell)
        label = f"length={nominal};depth={depth}x"
        for metric in ("runs_with_STRAND_BIASED_ALT_SUPPORT",
                       "runs_with_any_tier_affecting_flag", "accepted", "exact"):
            add(2, label, metric, cell[metric], len(rows), "test2_cells.tsv")
    _write_tsv(out_dir / "test2_cells.tsv", cells, list(cells[0]) if cells else [])
    if by_test[2]:
        add(2, "all", "runs_with_STRAND_BIASED_ALT_SUPPORT",
            sum(c["runs_with_STRAND_BIASED_ALT_SUPPORT"] for c in cells),
            sum(c["runs"] for c in cells), "test2_cells.tsv")

    # ---- test 3 ----------------------------------------------------------- #
    cells = []
    groups = OrderedDict()
    for r in by_test[3]:
        observed_shift = _has(r, "FRAME_LENGTH_SHIFT")
        observed_stops = _has(r, "INTERNAL_STOPS")
        true_shift = r.get("truth_length_mod3") != r.get("bait_length_mod3")
        true_stops = (r.get("truth_internal_stops_min_six_frames") or 0) > 0
        r["frame_flags_observed"] = ";".join(
            n for n, on in (("FRAME_LENGTH_SHIFT", observed_shift),
                            ("INTERNAL_STOPS", observed_stops)) if on)
        r["frame_flags_expected_from_truth"] = ";".join(
            n for n, on in (("FRAME_LENGTH_SHIFT", true_shift),
                            ("INTERNAL_STOPS", true_stops)) if on)
        r["frame_flags_agree_with_truth"] = (
            int(observed_shift == true_shift and observed_stops == true_stops)
            if int(r.get("sequence_returned") or 0) else "")
        # against the true reading frame (+1 of the bait): is a disrupted frame
        # alerted, and is an intact frame left without alert?
        r["frame_alert_matches_true_frame"] = (
            int((observed_shift or observed_stops) == (not r.get("truth_frame_intact")))
            if int(r.get("sequence_returned") or 0) else "")
        groups.setdefault((r["scenario"], r["depth"]), []).append(r)
    _write_tsv(out_dir / TABLES[3][0], by_test[3],
               TABLES[3][1] + ["frame_flags_observed", "frame_flags_expected_from_truth",
                               "frame_flags_agree_with_truth",
                               "frame_alert_matches_true_frame"])
    for (scenario, depth), rows in groups.items():
        returned = [r for r in rows if int(r.get("sequence_returned") or 0)]
        cell = {
            "scenario": scenario, "depth": depth, "runs": len(rows),
            "sequence_returned": len(returned),
            "exact": sum(int(r.get("exact") or 0) for r in rows),
            "accepted": sum(1 for r in rows if r.get("lr_result_disposition") == "PASS"),
            "accepted_not_exact": sum(1 for r in rows if r.get("lr_result_disposition") == "PASS"
                                      and not int(r.get("exact") or 0)),
            "sent_to_review": sum(1 for r in rows if r.get("lr_result_disposition") == "REVIEW"),
            "held": sum(1 for r in rows if r.get("lr_result_disposition") == "HOLD"),
            "tiers": _counts(r.get("lr_qc_confidence") for r in rows),
            "uncertain_bases_interior_median_range": _range(
                [_num(r.get("lr_interior_uncertain_base_count")) for r in returned]),
            "uncertain_bases_whole_median_range": _range(
                [_num(r.get("lr_uncertain_base_count")) for r in returned]),
            "true_frame_intact": _counts(r.get("truth_frame_intact") for r in rows),
            "frame_flags_expected_from_truth": _counts(
                r["frame_flags_expected_from_truth"] for r in rows),
            "frame_flags_observed": _counts(r["frame_flags_observed"] for r in returned),
            "frame_flags_agree_with_truth": sum(
                1 for r in returned if r["frame_flags_agree_with_truth"] == 1),
            "frame_alert_matches_true_frame": sum(
                1 for r in returned if r["frame_alert_matches_true_frame"] == 1),
            "not_exact_detail": "; ".join(
                f"{r['run_id']}:mm={r.get('cmp_mismatches')},ins={r.get('cmp_inserted_bp')},"
                f"del={r.get('cmp_deleted_bp')}" for r in returned if not int(r.get("exact") or 0)),
        }
        cells.append(cell)
        label = f"{scenario};depth={depth}x"
        for metric in ("exact", "accepted", "accepted_not_exact", "held",
                       "frame_flags_agree_with_truth"):
            add(3, label, metric, cell[metric],
                len(returned) if metric == "frame_flags_agree_with_truth" else len(rows),
                "test3_cells.tsv")
    _write_tsv(out_dir / "test3_cells.tsv", cells, list(cells[0]) if cells else [])

    # ---- tests 4 and 5: classification against the truth ------------------ #
    for r in by_test[4] + by_test[5]:
        truth = r.get("true_copies")
        count = _num(r.get("graph_context_count"))
        bound = _num(r.get("graph_context_lower_bound"))
        if r.get("cn_status") != "DONE":
            graph_class = "not_evaluated"
        elif count is not None:
            graph_class = "equals_true_copies" if count == truth else "exact_count_wrong"
        elif bound is not None:
            graph_class = ("lower_bound_correct" if bound <= truth else "lower_bound_wrong")
        else:
            graph_class = "no_count"
        r["graph_vs_truth"] = graph_class
        for label in ("cn", "cn_nograph"):
            kind = r.get(f"{label}_copy_number_kind")
            call = _num(r.get(f"{label}_copy_number_call"))
            lower = _num(r.get(f"{label}_copy_number_lower_bound"))
            if not kind:
                verdict = "not_evaluated"
            elif kind == "INTEGER_CONTEXT_COUNT":
                verdict = "integer_equals_truth" if call == truth else "integer_wrong"
            elif kind == "MEAN_DEPTH_DOSAGE":
                verdict = ("dosage_within_0.5_of_truth" if call is not None
                           and abs(call - truth) < 0.5 else "dosage_off_by_0.5_or_more")
            elif kind == "LOWER_BOUND":
                verdict = ("lower_bound_correct" if lower is not None and lower <= truth
                           else "lower_bound_wrong")
            else:
                verdict = "not_estimated"
            r[f"{label}_vs_truth"] = verdict
            r[f"{label}_supported_but_wrong"] = int(
                r.get(f"{label}_copy_number_status") == "SUPPORTED"
                and verdict in ("integer_wrong", "dosage_off_by_0.5_or_more"))
        dosage = _num(r.get("dr_dosage_estimate"))
        r["dosage_abs_error"] = round(abs(dosage - truth), 3) if dosage is not None else ""
        r["depth_call_agrees_with_truth"] = (
            int((r.get("dr_depth_call") == "MULTICOPY_DEPTH") == (truth >= 2))
            if r.get("dr_depth_call") else "")
    extra = ["graph_vs_truth", "cn_vs_truth", "cn_supported_but_wrong",
             "cn_nograph_vs_truth", "cn_nograph_supported_but_wrong", "dosage_abs_error",
             "depth_call_agrees_with_truth"]
    for test in (4, 5):
        _write_tsv(out_dir / TABLES[test][0], by_test[test], TABLES[test][1] + extra)

    def consensus_text(rows, label):
        return _counts(
            f"{r.get(label + '_copy_number_call')}|{r.get(label + '_copy_number_kind')}|"
            f"{r.get(label + '_copy_number_method')}|{r.get(label + '_copy_number_status')}"
            for r in rows)

    cells = []
    groups = OrderedDict()
    for r in by_test[4]:
        groups.setdefault((r.get("background", "plain"), r["scenario"]), []).append(r)
    for (background, scenario), rows in groups.items():
        cell = {
            "background": background,
            "scenario": scenario, "true_copies": rows[0]["true_copies"], "runs": len(rows),
            "truth": rows[0].get("truth", ""),
            "discovery_intervals": _counts(r.get("cn_discovery_intervals") for r in rows),
            "left_context_count": _counts(r.get("left_context_count") for r in rows),
            "right_context_count": _counts(r.get("right_context_count") for r in rows),
            "graph_context_count": _counts(r.get("graph_context_count") for r in rows),
            "graph_context_lower_bound": _counts(r.get("graph_context_lower_bound") for r in rows),
            "graph_context_status": _counts(r.get("graph_context_status") for r in rows),
            "graph_context_flags": _counts(r.get("graph_context_flags") for r in rows),
            "graph_vs_truth": _counts(r.get("graph_vs_truth") for r in rows),
            "depth_call": _counts(r.get("dr_depth_call") for r in rows),
            "depth_call_agrees_with_truth": sum(
                1 for r in rows if r.get("depth_call_agrees_with_truth") == 1),
            "ratio_median_range": _range([_num(r.get("dr_ratio")) for r in rows]),
            "dosage_median_range": _range([_num(r.get("dr_dosage_estimate")) for r in rows]),
            "dosage_status": _counts(r.get("dr_dosage_status") for r in rows),
            "ambiguity_index_median_range": _range(
                [_num(r.get("dr_ambiguity_index")) for r in rows]),
            "consensus_call_kind_method_status": consensus_text(rows, "cn"),
            "consensus_vs_truth": _counts(r.get("cn_vs_truth") for r in rows),
            "consensus_supported_but_wrong": sum(r.get("cn_supported_but_wrong", 0) for r in rows),
            "locus_recon_dispositions": _counts(r.get("lr_result_disposition") for r in rows),
            "locus_recon_exact": sum(int(r.get("exact") or 0) for r in rows),
        }
        cells.append(cell)
        for metric in ("graph_vs_truth", "depth_call", "consensus_vs_truth",
                       "consensus_supported_but_wrong", "dosage_median_range"):
            add(4, f"background={background};{scenario}", metric, cell[metric], len(rows),
                "test4_cells.tsv")
    _write_tsv(out_dir / "test4_cells.tsv", cells, list(cells[0]) if cells else [])

    cells = []
    groups = OrderedDict()
    for r in by_test[5]:
        groups.setdefault((r.get("arm", ""), r["true_copies"], r["overlap_pct"]), []).append(r)
    for (arm, copies, overlap), rows in groups.items():
        cell = {
            "arm": arm,
            "true_copies": copies, "overlap_pct": overlap, "runs": len(rows),
            "cut_status": _counts(r.get("cut_status") for r in rows),
            "discovery_intervals": _counts(r.get("cn_discovery_intervals") for r in rows),
            "gene_overlap_fraction": _counts(r.get("dr_gene_overlap_fraction") for r in rows),
            "region_ratios": " | ".join(str(r.get("dr_region_ratios", "")) for r in rows),
            "dosage_median_range": _range([_num(r.get("dr_dosage_estimate")) for r in rows]),
            "dosage_ci_low_min": min([v for v in (_num(r.get("dr_dosage_ci_low")) for r in rows)
                                      if v is not None], default=""),
            "dosage_ci_high_max": max([v for v in (_num(r.get("dr_dosage_ci_high")) for r in rows)
                                       if v is not None], default=""),
            "geometry_expected_dosage": _counts(r.get("dr_geometry_expected_dosage") for r in rows),
            "dosage_status": _counts(r.get("dr_dosage_status") for r in rows),
            "dosage_flags": _counts(r.get("dr_dosage_flags") for r in rows),
            "ambiguity_index_median_range": _range(
                [_num(r.get("dr_ambiguity_index")) for r in rows]),
            "depth_call": _counts(r.get("dr_depth_call") for r in rows),
            "consensus_without_graph": consensus_text(rows, "cn_nograph"),
            "consensus_status_without_graph": _counts(
                r.get("cn_nograph_copy_number_status") for r in rows),
            "consensus_with_graph": consensus_text(rows, "cn"),
            "consensus_status_with_graph": _counts(r.get("cn_copy_number_status") for r in rows),
            "graph_context_count": _counts(r.get("graph_context_count") for r in rows),
            "supported_but_wrong_without_graph": sum(
                r.get("cn_nograph_supported_but_wrong", 0) for r in rows),
            "supported_but_wrong_with_graph": sum(r.get("cn_supported_but_wrong", 0) for r in rows),
            "locus_recon_dispositions": _counts(r.get("lr_result_disposition") for r in rows),
            "locus_recon_exact": sum(int(r.get("exact") or 0) for r in rows),
        }
        cells.append(cell)
        label = f"arm={arm};overlap={overlap}%"
        for metric in ("dosage_median_range", "geometry_expected_dosage", "dosage_status",
                       "consensus_status_without_graph", "consensus_status_with_graph"):
            add(5, label, metric, cell[metric], len(rows), "test5_cells.tsv")
    _write_tsv(out_dir / "test5_cells.tsv", cells, list(cells[0]) if cells else [])

    separation = []
    two_copy_arms = sorted({r.get("arm", "") for r in by_test[5] if r["true_copies"] == 2})
    for two_arm, overlap in ((arm, ov) for arm in two_copy_arms for ov in T5_OVERLAPS):
        one_rows = [r for r in by_test[5]
                    if r["true_copies"] == 1 and r["overlap_pct"] == overlap]
        two_rows = [r for r in by_test[5]
                    if r.get("arm", "") == two_arm and r["overlap_pct"] == overlap]
        one = [_num(r.get("dr_dosage_estimate")) for r in one_rows]
        two = [_num(r.get("dr_dosage_estimate")) for r in two_rows]
        one = [v for v in one if v is not None]
        two = [v for v in two if v is not None]
        if not one or not two:
            continue
        one_ci = [_num(r.get("dr_dosage_ci_high")) for r in one_rows]
        two_ci = [_num(r.get("dr_dosage_ci_low")) for r in two_rows]
        one_ci = [v for v in one_ci if v is not None]
        two_ci = [v for v in two_ci if v is not None]
        row = {
            "one_copy_arm": "c1_identical_cut", "two_copy_arm": two_arm,
            "overlap_pct": overlap,
            "one_copy_runs": len(one), "two_copy_runs": len(two),
            "one_copy_dosage_min": min(one), "one_copy_dosage_max": max(one),
            "two_copy_dosage_min": min(two), "two_copy_dosage_max": max(two),
            "dosage_ranges_disjoint": int(max(one) < min(two)),
            "gap_between_ranges": round(min(two) - max(one), 3),
            "one_copy_ci_high_max": max(one_ci) if one_ci else "",
            "two_copy_ci_low_min": min(two_ci) if two_ci else "",
            "intervals_disjoint": (int(max(one_ci) < min(two_ci)) if one_ci and two_ci else ""),
            "ratio_two_over_one_medians": round(
                statistics.median(two) / statistics.median(one), 3)
            if statistics.median(one) > 0 else "",
        }
        separation.append(row)
        add(5, f"c1_identical_cut vs {two_arm};overlap={overlap}%",
            "dosage_ranges_disjoint_one_vs_two_copies",
            row["dosage_ranges_disjoint"], "", "test5_separation.tsv")
    _write_tsv(out_dir / "test5_separation.tsv", separation,
               list(separation[0]) if separation else [])

    # ---- block sizes ------------------------------------------------------ #
    block_rows = []
    for r in by_test[4] + by_test[5]:
        for row in r.get("block_rows", []) or []:
            block_rows.append({"run_id": r["run_id"], "test": r["test"],
                               "background": r.get("background", ""),
                               "scenario": r["scenario"], **row})
    block_columns = ["run_id", "test", "background", "scenario", "block_bp", "true_copies",
                     "assembly_region_count", "ratio", "ratio_ci_low", "ratio_ci_high",
                     "ratio_ci_width", "ratio_ci_contains_true_copies", "depth_call",
                     "dosage_estimate", "dosage_ci_low", "dosage_ci_high", "dosage_ci_width",
                     "dosage_ci_contains_true_copies", "dosage_status",
                     "bootstrap_locus_blocks", "bootstrap_backbone_blocks_total",
                     "bootstrap_backbone_blocks_used"]
    _write_tsv(out_dir / "block_size_tests4_5.tsv", block_rows, block_columns)
    cells = []
    groups = OrderedDict()
    for row in block_rows:
        groups.setdefault(
            (row["test"], row["background"], row["scenario"], row["block_bp"]), []).append(row)
    for (test, background, scenario, block), rows in groups.items():
        widths = [_num(r["dosage_ci_width"]) for r in rows]
        contains = [r["dosage_ci_contains_true_copies"] for r in rows
                    if r["dosage_ci_contains_true_copies"] != ""]
        cells.append({
            "test": test, "background": background, "scenario": scenario,
            "block_bp": block, "runs": len(rows),
            "true_copies": rows[0]["true_copies"],
            "dosage_ci_width_median_range": _range(widths),
            "dosage_ci_contains_true_copies": sum(contains),
            "ratio_ci_width_median_range": _range([_num(r["ratio_ci_width"]) for r in rows]),
            "ratio_ci_contains_true_copies": sum(r["ratio_ci_contains_true_copies"] for r in rows),
            "depth_call": _counts(r["depth_call"] for r in rows),
        })
    _write_tsv(out_dir / "block_size_cells.tsv", cells, list(cells[0]) if cells else [])
    overall = OrderedDict()
    for row in block_rows:
        overall.setdefault((f"{row['test']}:{row['background']}", row["block_bp"]), []).append(row)
    block_summary = []
    for (test, block), rows in overall.items():
        contains = [r["dosage_ci_contains_true_copies"] for r in rows
                    if r["dosage_ci_contains_true_copies"] != ""]
        calls_default = {r["run_id"]: r["depth_call"] for r in block_rows
                         if r["block_bp"] == 500}
        changed = sum(1 for r in rows if calls_default.get(r["run_id"]) != r["depth_call"])
        entry = {
            "test": test, "block_bp": block, "runs": len(rows),
            "dosage_ci_width_median_range": _range([_num(r["dosage_ci_width"]) for r in rows]),
            "dosage_ci_contains_true_copies": sum(contains),
            "dosage_ci_evaluated": len(contains),
            "ratio_ci_width_median_range": _range([_num(r["ratio_ci_width"]) for r in rows]),
            "depth_call_differs_from_500bp": changed,
        }
        block_summary.append(entry)
        add(f"{test}-blocks", f"block={block}bp", "dosage_ci_contains_true_copies",
            entry["dosage_ci_contains_true_copies"], entry["dosage_ci_evaluated"],
            "block_size_summary.tsv")
        add(f"{test}-blocks", f"block={block}bp", "dosage_ci_width_median_range",
            entry["dosage_ci_width_median_range"], len(rows), "block_size_summary.tsv")
        add(f"{test}-blocks", f"block={block}bp", "depth_call_differs_from_500bp",
            changed, len(rows), "block_size_summary.tsv")
    _write_tsv(out_dir / "block_size_summary.tsv", block_summary,
               list(block_summary[0]) if block_summary else [])

    # ---- marginal totals, counted from the per-run tables ------------------ #
    def tally(test, label, rows, source, extra=()):
        n = len(rows)
        if not n:
            return
        accepted = [r for r in rows if r.get("lr_result_disposition") == "PASS"]
        held = [r for r in rows if r.get("lr_result_disposition") == "HOLD"]
        values = [
            ("sequence_returned", sum(int(r.get("sequence_returned") or 0) for r in rows)),
            ("exact", sum(int(r.get("exact") or 0) for r in rows)),
            ("accepted", len(accepted)),
            ("accepted_not_exact", sum(1 for r in accepted if not int(r.get("exact") or 0))),
            ("sent_to_review",
             sum(1 for r in rows if r.get("lr_result_disposition") == "REVIEW")),
            ("held", len(held)),
            ("held_and_exact", sum(1 for r in held if int(r.get("exact") or 0))),
            ("runs_with_STRAND_BIASED_ALT_SUPPORT",
             sum(1 for r in rows if _has(r, "STRAND_BIASED_ALT_SUPPORT"))),
            ("runs_with_MIXED_SITE_BELOW_MIXTURE_CALL",
             sum(1 for r in rows if _has(r, "MIXED_SITE_BELOW_MIXTURE_CALL"))),
        ]
        for metric, value in values + [(m, f(rows)) for m, f in extra]:
            add(test, label, metric, value, n, source)

    source = TABLES[1][0]
    mixture_extra = (
        ("held_as_mixture", lambda rows: sum(
            1 for r in rows if r.get("lr_result_disposition") == "HOLD"
            and r.get("lr_mixture_detected") == "true")),
        ("accepted_with_MIXED_SITE_BELOW_MIXTURE_CALL", lambda rows: sum(
            1 for r in rows if r.get("lr_result_disposition") == "PASS"
            and _has(r, "MIXED_SITE_BELOW_MIXTURE_CALL"))),
        ("accepted_without_any_mixed_site_flag", lambda rows: sum(
            1 for r in rows if r.get("lr_result_disposition") == "PASS"
            and not _has(r, "MIXED_SITE_BELOW_MIXTURE_CALL"))),
        ("exact_major_allele", lambda rows: sum(
            1 for r in rows if r.get("matches_allele") in ("major", "allele"))),
        ("exact_minor_allele", lambda rows: sum(
            1 for r in rows if r.get("matches_allele") == "minor")),
    )
    mixtures = [r for r in by_test[1] if r["scenario"] == "mixture"]
    for sites in T1_SITES:
        tally("1-total", f"sites={sites};all fractions",
              [r for r in mixtures if r.get("sites") == sites], source, mixture_extra)
    for fraction in T1_FRACTIONS:
        tally("1-total", f"minor={fraction}%;all site numbers",
              [r for r in mixtures if r.get("minor_fraction_pct") == fraction], source,
              mixture_extra)
    tally("1-total", "all mixtures", mixtures, source, mixture_extra)
    tally("1-total", "pure controls",
          [r for r in by_test[1] if r["scenario"] == "pure_control"], source, mixture_extra)

    source = TABLES[2][0]
    tier_flag = (("runs_with_any_tier_affecting_flag", lambda rows: sum(
        1 for r in rows if any(n and not n.startswith(DESCRIPTIVE)
                               for n in (r.get("lr_flag_names") or "").split(";")))),)
    for depth in T2_DEPTHS:
        tally("2-total", f"depth={depth}x;all lengths",
              [r for r in by_test[2] if r["depth"] == depth], source, tier_flag)
    for length in T2_LENGTHS:
        tally("2-total", f"length={T2_NOMINAL[length]};all depths",
              [r for r in by_test[2] if r["locus_bp"] == length], source, tier_flag)
    tally("2-total", "all", by_test[2], source, tier_flag)

    source = TABLES[3][0]
    frame_extra = (
        ("frame_flags_agree_with_truth", lambda rows: sum(
            1 for r in rows if r.get("frame_flags_agree_with_truth") == 1)),
        ("frame_alert_matches_true_frame", lambda rows: sum(
            1 for r in rows if r.get("frame_alert_matches_true_frame") == 1)),
        ("true_frame_disrupted", lambda rows: sum(
            1 for r in rows if r.get("truth_frame_intact") == 0)),
        ("true_frame_disrupted_and_held", lambda rows: sum(
            1 for r in rows if r.get("truth_frame_intact") == 0
            and r.get("lr_result_disposition") == "HOLD")),
    )
    for scenario in T3_SCENARIOS:
        tally("3-total", f"{scenario};all depths",
              [r for r in by_test[3] if r["scenario"] == scenario], source, frame_extra)
    for depth in T3_DEPTHS:
        tally("3-total", f"depth={depth}x;all scenarios",
              [r for r in by_test[3] if r["depth"] == depth], source, frame_extra)
    tally("3-total", "all", by_test[3], source, frame_extra)

    source = TABLES[4][0]
    copy_extra = (
        ("graph_vs_truth", lambda rows: _counts(r.get("graph_vs_truth") for r in rows)),
        ("graph_context_status",
         lambda rows: _counts(r.get("graph_context_status") for r in rows)),
        ("depth_call", lambda rows: _counts(r.get("dr_depth_call") for r in rows)),
        ("depth_call_agrees_with_truth", lambda rows: sum(
            1 for r in rows if r.get("depth_call_agrees_with_truth") == 1)),
        ("dosage_status", lambda rows: _counts(r.get("dr_dosage_status") for r in rows)),
        ("dosage_within_0.5_of_truth", lambda rows: sum(
            1 for r in rows if r.get("dosage_abs_error") != ""
            and r.get("dosage_abs_error") < 0.5)),
        ("consensus_status",
         lambda rows: _counts(r.get("cn_copy_number_status") for r in rows)),
        ("consensus_vs_truth", lambda rows: _counts(r.get("cn_vs_truth") for r in rows)),
        ("consensus_supported_but_wrong",
         lambda rows: sum(r.get("cn_supported_but_wrong", 0) for r in rows)),
        ("consensus_status_without_graph",
         lambda rows: _counts(r.get("cn_nograph_copy_number_status") for r in rows)),
        ("consensus_without_graph_supported_but_wrong",
         lambda rows: sum(r.get("cn_nograph_supported_but_wrong", 0) for r in rows)),
    )
    for background in T4_BACKGROUNDS:
        arm = [r for r in by_test[4] if r.get("background") == background]
        tally("4-total", f"background={background};one copy",
              [r for r in arm if r["true_copies"] == 1], source, copy_extra)
        tally("4-total", f"background={background};two to four identical copies",
              [r for r in arm if r["true_copies"] > 1 and not r["divergent"]], source,
              copy_extra)
        tally("4-total", f"background={background};two to four copies, 1% divergent",
              [r for r in arm if r["true_copies"] > 1 and r["divergent"]], source, copy_extra)
        tally("4-total", f"background={background};all", arm, source, copy_extra)
    tally("4-total", "all", by_test[4], source, copy_extra)

    source = TABLES[5][0]
    for arm_name in sorted({r.get("arm", "") for r in by_test[5]}):
        arm = [r for r in by_test[5] if r.get("arm", "") == arm_name]
        tally("5-total", f"arm={arm_name};all overlaps", arm, source, copy_extra)
        tally("5-total", f"arm={arm_name};overlap>=45%",
              [r for r in arm if r["overlap_pct"] >= 45], source, copy_extra)
    tally("5-total", "all", by_test[5], source, copy_extra)

    status = Counter((r["test"], r.get("run_status")) for r in records)
    for (test, state), count in sorted(status.items(), key=lambda item: str(item)):
        add(test, "all", f"runs_with_status_{state}", count,
            sum(1 for r in records if r["test"] == test), "targeted_manifest.tsv")
    _write_tsv(out_dir / "targeted_summary.tsv", summary,
               ["test", "cell", "metric", "value", "denominator", "source_file"])


def collect_sequences(records: list[dict], work_dir: Path, out_dir: Path) -> None:
    """Truth alleles, baits and returned sequences of every run, per test."""
    for test in (1, 2, 3, 4, 5):
        truth_records, bait_records, returned_records = [], [], []
        for r in records:
            if r["test"] != test:
                continue
            base = work_dir / "runs" / f"test{test}" / r["run_id"]
            if (base / "truth" / "truth_alleles.fasta").exists():
                truth_records += read_fasta(base / "truth" / "truth_alleles.fasta")
            if (base / "truth" / "bait.fasta").exists():
                bait_records += [(f"{r['run_id']}|bait", s)
                                 for _n, s in read_fasta(base / "truth" / "bait.fasta")]
            if (base / "returned.fasta").exists():
                returned_records += read_fasta(base / "returned.fasta")
        seq_dir = out_dir / "sequences"
        seq_dir.mkdir(exist_ok=True)
        write_fasta(seq_dir / f"test{test}_truth_alleles.fasta", truth_records)
        write_fasta(seq_dir / f"test{test}_baits.fasta", bait_records)
        write_fasta(seq_dir / f"test{test}_returned.fasta", returned_records)


def provenance(args, records: list[dict], out_dir: Path, wall_seconds: float) -> None:
    import Bio
    import locus_recon
    package_root = Path(locus_recon.__file__).resolve().parent.parent
    digest = source_sha256(package_root)
    tools = {
        "python": [sys.executable, "--version"],
        "locus-recon": ["locus-recon", "--version"],
        "locus-recon-depth-ratio": ["locus-recon-depth-ratio", "--version"],
        "locus-recon-copy-number": ["locus-recon-copy-number", "--version"],
        "spades": ["spades.py", "--version"],
        "bwa-mem2": ["bwa-mem2", "version"],
        "samtools": ["samtools", "--version"],
        "blastn": ["blastn", "-version"],
        "makeblastdb": ["makeblastdb", "-version"],
    }
    record = {
        "schema_version": "1.0",
        "workflow": "targeted-tests-simulated-reads",
        "protocol": "validation/external-evaluation/PROTOCOL.md, section 3",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "software": {
            "name": "locus-recon", "version": locus_recon.VERSION,
            "package_sha256": digest, "expected_package_sha256": EXPECTED_PACKAGE_DIGEST,
            "package_sha256_matches": digest == EXPECTED_PACKAGE_DIGEST,
            "package_location": "installed package (site-packages or checkout)",
        },
        "runner": {
            "path": "validation/run_targeted_tests.py",
            "sha256": sha256_file(Path(__file__).resolve()),
            "generator_sha256": sha256_file(VALIDATION_DIR / "generate_mock_dataset.py"),
            "common_sha256": sha256_file(VALIDATION_DIR / "common.py"),
        },
        "tools": {name: command_version(cmd)["version"] for name, cmd in tools.items()},
        "libraries": {"numpy": np.__version__, "biopython": Bio.__version__,
                      "zlib": zlib.ZLIB_VERSION},
        "platform": platform.platform(),
        "seeds": {
            "base_seed": BASE_SEED,
            "rule": "seed = int.from_bytes(sha256(f'{BASE_SEED}|{run_id}|{purpose}')[:8], 'big') >> 1; "
                    "purpose is 'genome' or 'reads'; numpy.random.default_rng(seed)",
            "bootstrap_seed": BOOTSTRAP_SEED, "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
            "per_run_seeds": "targeted_manifest.tsv (seed_genome, seed_reads)",
        },
        "simulation": {
            "read_length": READ_LENGTH, "fragment_mean": INSERT_SIZE, "fragment_sd": INSERT_SD,
            "fragment_min": INSERT_MIN, "fragment_max": INSERT_MAX,
            "substitution_error_rate": ERROR_RATE, "base_quality": "I (Phred 40)",
            "fragment_starts": "uniform", "orientation": "random", "background_bp": BACKGROUND_BP,
        },
        "commands": {
            "draft": "spades.py --isolate -1 R1 -2 R2 -o DIR --threads T --memory M",
            "locus_recon": "locus-recon --samplesheet S --main-output-dir DIR --bait BAIT "
                           f"--locus {LOCUS} --threads T --memory-per-sample M",
            "discovery_for_copy_number": "makeblastdb -in DRAFT -dbtype nucl; blastn -query BAIT "
                                         "-db DRAFT -outfmt '6 qseqid sseqid pident length qstart "
                                         "qend sstart send bitscore' -max_target_seqs 20 "
                                         "-perc_identity 80",
            "depth_ratio": "locus-recon-depth-ratio --bam SAMPLE/full_map.sorted.bam "
                           "{--locus REGION --gene-span A-B}... --gene-length L "
                           "--reference DRAFT --tsv OUT",
            "copy_number": "locus-recon-copy-number (same depth options) --gfa "
                           "DRAFT_SPADES/assembly_graph_after_simplification.gfa --bait BAIT "
                           "--output-json OUT.json --output-tsv OUT.tsv",
            "copy_number_without_graph": "locus-recon-copy-number (same depth options) "
                                         "--output-json OUT.json --output-tsv OUT.tsv",
            "block_sizes": "locus_recon.depth_ratio.estimate_locus_copy_number(bam, regions, "
                           "gene_spans=spans, gene_length=L, block=B, n_boot=1000, "
                           "seed=20260817) for B in (100, 250, 500, 1000)",
        },
        "parameters": {
            "jobs": args.jobs, "threads_per_run": args.threads_per_run,
            "memory_gb_per_run": args.memory_gb, "replicate_offset": args.replicate_offset,
            "replicate_cap": args.replicate_cap, "tests": args.tests,
        },
        "runs": {"planned": len(records),
                 "by_status": dict(Counter(r.get("run_status") for r in records))},
        "wall_seconds_this_invocation": round(wall_seconds, 1),
        "input_checksums": "input_checksums.tsv (SHA-256 of reads, bait, truth, draft, graph "
                           "and alignment of every run)",
    }
    write_json(out_dir / "PROVENANCE.json", record)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work-dir", type=Path, required=True,
                        help="heavy per-run outputs (reads, drafts, alignments)")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="tables, sequences and provenance")
    parser.add_argument("--jobs", type=int, default=8, help="runs in parallel")
    parser.add_argument("--threads-per-run", type=int, default=1)
    parser.add_argument("--memory-gb", type=int, default=8)
    parser.add_argument("--tests", default="1,2,3,4,5")
    parser.add_argument("--only", default="", help="comma-separated run identifiers")
    parser.add_argument("--replicate-offset", type=int, default=0,
                        help="pilot runs only: shifts the replicate numbers (and the seeds)")
    parser.add_argument("--replicate-cap", type=int, default=None,
                        help="pilot runs only: replicates per cell")
    parser.add_argument("--aggregate-only", action="store_true")
    parser.add_argument("--allow-digest-mismatch", action="store_true")
    args = parser.parse_args()

    import locus_recon
    digest = source_sha256(Path(locus_recon.__file__).resolve().parent.parent)
    if (digest != EXPECTED_PACKAGE_DIGEST or locus_recon.VERSION != EXPECTED_VERSION) \
            and not args.allow_digest_mismatch:
        raise SystemExit(
            f"installed locus-recon {locus_recon.VERSION} has package digest {digest}; "
            f"the protocol freezes {EXPECTED_VERSION} / {EXPECTED_PACKAGE_DIGEST}")
    for tool in ("spades.py", "locus-recon", "locus-recon-depth-ratio",
                 "locus-recon-copy-number", "blastn", "makeblastdb", "samtools"):
        if not shutil.which(tool):
            raise FileNotFoundError(f"{tool} is required on PATH")

    tests = {int(t) for t in args.tests.split(",") if t}
    runs = [r for r in build_manifest(args.replicate_offset, args.replicate_cap)
            if r["test"] in tests]
    if args.only:
        wanted = set(args.only.split(","))
        runs = [r for r in runs if r["run_id"] in wanted]
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work_dir = args.work_dir.resolve()
    started = time.time()
    records: dict[str, dict] = {}
    if args.aggregate_only:
        for run in runs:
            path = work_dir / "runs" / f"test{run['test']}" / run["run_id"] / "result.json"
            if path.exists():
                records[run["run_id"]] = json.loads(path.read_text())
            else:
                records[run["run_id"]] = {**run, "run_status": "NOT_RUN"}
    else:
        # heavier runs first, so the pool stays busy at the end
        order = sorted(runs, key=lambda r: -(r["depth"] * (2 if r["test"] in (4, 5) else 1)))
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(process_run, run, str(work_dir), args.threads_per_run,
                                   args.memory_gb): run for run in order}
            for count, future in enumerate(as_completed(futures), start=1):
                run = futures[future]
                try:
                    record = future.result()
                except Exception as exc:
                    record = {**run, "run_status": "RUNNER_ERROR",
                              "runner_error": f"{type(exc).__name__}: {exc}"}
                records[run["run_id"]] = record
                print(f"[{count}/{len(order)}] {run['run_id']} {record.get('run_status')} "
                      f"{record.get('lr_result_disposition', '')} exact={record.get('exact', '')} "
                      f"{record.get('elapsed_sec', '')}s {record.get('runner_error', '')}",
                      flush=True)
    ordered = [records[r["run_id"]] for r in runs]
    aggregate(ordered, args.out_dir)
    collect_sequences(ordered, work_dir, args.out_dir)
    provenance(args, ordered, args.out_dir, time.time() - started)
    print(f"tables written to {args.out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
