"""
Shared constants, exceptions, logging setup, dependency checking, and
command execution helpers for Locus-Recon.
"""

import logging
import re
import shlex
import shutil
import subprocess
import sys
from typing import Tuple

from . import PROGRAM, VERSION  # noqa: F401 – re-exported for convenience

# ==============================================================================
# Constants
# ==============================================================================

BLAST_FMT_STEP1 = (
    "6 qseqid sseqid qstart qend sstart send length pident qcovhsp bitscore evalue"
)
BLAST_FMT_STEP5 = (
    "6 qseqid sseqid qstart qend sstart send length pident qcovhsp "
    "bitscore evalue qlen slen"
)

# Stop codons (standard bacterial genetic code 11)
STOP_CODONS = {"TAA", "TAG", "TGA"}

# QC confidence thresholds
QC_THRESHOLDS = {
    "high": {
        "length_iqr_factor": 1.0,
        "min_identity": 97.0,
        "min_qcov": 95.0,
        "max_n_fraction": 0.001,
        "max_gc_deviation": 3.0,
    },
    "medium": {
        "length_iqr_factor": 2.0,
        "min_identity": 93.0,
        "min_qcov": 90.0,
        "max_n_fraction": 0.005,
        "max_gc_deviation": 5.0,
    },
    "low": {
        "length_iqr_factor": 3.0,
        "min_identity": 85.0,
        "min_qcov": 80.0,
        "max_n_fraction": 0.01,
        "max_gc_deviation": 8.0,
    },
    "novel_override_identity": 97.0,
}

REMAP_THRESHOLDS = {
    "high":   {"min_mean_depth": 20.0, "min_breadth": 99.0, "max_pct_bases_lt5": 1.0},
    "medium": {"min_mean_depth": 10.0, "min_breadth": 97.0, "max_pct_bases_lt5": 5.0},
    "low":    {"min_mean_depth":  5.0, "min_breadth": 90.0, "max_pct_bases_lt5": 15.0},
}

# Per-base certainty model.
BASE_CERTAINTY_THRESHOLDS = {
    # A position is descriptively "low coverage" below this filtered depth.
    # This is a coverage statement only and does not imply base uncertainty.
    "low_coverage_depth": 5,
    # A single filtered observation cannot corroborate itself.
    "min_filtered_depth": 2,
    # Below this filtered depth, strand-restricted evidence is not sufficient.
    "bidirectional_required_below_depth": 5,
    # Alternative-allele fraction at which the reported base is not resolved.
    "max_alt_fraction": 0.20,
}
# Tier gates on interior_uncertain_base_fraction (fraction, not percent).
UNCERTAINTY_THRESHOLDS = {
    "high":   0.01,
    "medium": 0.05,
    "low":    0.15,
}
# Flags that describe a run but do not, on their own, constrain the tier.
DESCRIPTIVE_FLAG_PREFIXES = (
    "PATCHY_REMAP_SUPPORT",
    # A length that is not a multiple of three in a reconstruction already
    # known to be clipped by a contig boundary.  The truncation is scored by
    # its own flag; this one only records that the frame arithmetic follows
    # from it.
    "FRAME_LENGTH_SHIFT_TRUNCATED",
    # Internal stop codons that a scaffold placeholder's frameshift accounts
    # for.  INTERNAL_STOPS itself is still raised, because the delivered
    # sequence does contain them; this one records where they come from.
    "PLACEHOLDER_FRAMESHIFT_EXPLAINS_STOPS",
    # Graph evidence about a scaffold placeholder the local assembler inserted.
    # Reported so the reason for a held call is legible; it does not move the
    # tier, which remains a statement about read support for reported bases.
    "PLACEHOLDER_JUNCTION_",
    "MODEST_REMAP_PATCHINESS",
    "REMAP_PATCHINESS_BELOW_HIGH",
    # Catalogue relationship: how far the reconstruction sits from the nearest
    # curated allele.  Reported for interpretation and never used to constrain
    # sequence confidence, which is a statement about read support.
    "CATALOGUE_",
    "POSSIBLE_NOVEL_ALLELE",
)
# A projected locus span clipped by a contig boundary loses real query
# sequence, but terminal alignments routinely stop a base or two short of the
# query end for reasons that carry no completeness information.  Only clips at
# or above this size are treated as evidence that the allele is truncated.
SPAN_CLIP_THRESHOLDS = {
    "min_flag_bp": 10,
}
PER_BASE_THRESHOLDS = {
    "high": {
        "min_mean_base_quality": 30.0,
        "min_mean_mapping_quality": 40.0,
        "min_strand_balance": 40.0,
    },
    "medium": {
        "min_mean_base_quality": 25.0,
        "min_mean_mapping_quality": 30.0,
        "min_strand_balance": 25.0,
    },
    "low": {
        "min_mean_base_quality": 20.0,
        "min_mean_mapping_quality": 20.0,
        "min_strand_balance": 10.0,
    },
}

# ==============================================================================
# Exception
# ==============================================================================


class PipelineStepError(Exception):
    """Raised when a pipeline step fails for a single sample (non-fatal to batch)."""
    pass


class DependencyError(RuntimeError):
    """Raised when a required external executable is unavailable or too old."""

    pass


# ==============================================================================
# Module-level logger
# ==============================================================================

log = logging.getLogger(PROGRAM)

_LOG_FORMATTER = logging.Formatter(
    "[%(asctime)s] [%(levelname)-7s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def setup_logging(log_file: str, verbose: bool) -> None:
    """Configure the root logger with a file handler (DEBUG+) and a console
    handler (DEBUG if verbose, else WARNING+). Called once by the main process.

    Args:
        log_file: Path to the batch log file.
        verbose:  If True, console also shows DEBUG messages.
    """
    log.setLevel(logging.DEBUG)
    for handler in log.handlers:
        handler.close()
    log.handlers.clear()

    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(_LOG_FORMATTER)
    log.addHandler(fh)

    sh = logging.StreamHandler(sys.stderr)
    sh.setLevel(logging.DEBUG if verbose else logging.WARNING)
    sh.setFormatter(_LOG_FORMATTER)
    log.addHandler(sh)

    log.debug(f"Log file initialized: {log_file}")
    log.debug(f"Verbose mode: {'ON' if verbose else 'OFF'}")


def configure_sample_logger(sample_log_path: str, mode: str = "a") -> None:
    """Replace all handlers on the module logger with a single per-sample
    file handler. Used in parallel worker processes to avoid log contention.

    Args:
        sample_log_path: Path to the sample-specific log file.
    """
    logger = logging.getLogger(PROGRAM)
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(sample_log_path, mode=mode)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(_LOG_FORMATTER)
    logger.addHandler(fh)


def add_sample_log_handler(sample_log_path: str) -> logging.Handler:
    """Attach a per-sample file handler while retaining the batch handlers."""
    handler = logging.FileHandler(sample_log_path, mode="a")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(_LOG_FORMATTER)
    log.addHandler(handler)
    return handler


# ==============================================================================
# Dependency checking
# ==============================================================================


def check_dependencies(use_progress: bool) -> Tuple[dict, str]:
    """Verify all required CLI tools are available in $PATH.

    Args:
        use_progress: Whether tqdm progress bars are requested.

    Returns:
        tools:        Mapping of tool names to absolute paths.
        aligner_name: Either 'bwa-mem2' or 'bwa'.
    """
    log.info("Checking external dependencies...")

    tools: dict = {
        "blastn": None, "makeblastdb": None,
        "samtools": None, "spades.py": None,
    }
    if use_progress:
        tools["tqdm"] = None

    # Prefer bwa-mem2, fall back to bwa
    aligner_path = shutil.which("bwa-mem2")
    aligner_name = "bwa-mem2"
    if not aligner_path:
        aligner_path = shutil.which("bwa")
        aligner_name = "bwa"
    if not aligner_path:
        raise DependencyError(
            "No supported aligner found (bwa-mem2 or bwa). "
            "Create the environment with: mamba env create -f environment.yml"
        )

    tools[aligner_name] = aligner_path
    log.info(f"  Aligner : {aligner_name} -> {aligner_path}")

    missing = []
    for name in list(tools.keys()):
        if tools[name] is not None:
            continue
        path = shutil.which(name)
        if not path and name == "spades.py":
            path = shutil.which("spades")
        if path:
            tools[name] = path
            log.info(f"  Tool    : {name} -> {path}")
        else:
            if name == "tqdm":
                log.warning(f"  Optional: {name} not found -- progress bars disabled.")
                del tools[name]
            else:
                missing.append(name)

    if missing:
        raise DependencyError(
            "Missing required tools: " + ", ".join(missing) + ". "
            "Create the environment with: mamba env create -f environment.yml"
        )

    _check_samtools_version(tools["samtools"])
    log.info("All dependencies satisfied.")
    return tools, aligner_name


def _check_samtools_version(samtools_path: str) -> None:
    """Verify samtools >= 1.12 (required for 'samtools view -N' in mate rescue).

    Args:
        samtools_path: Absolute path to the samtools executable.
    """
    try:
        version_text = ""
        match = None
        for command in ([samtools_path, "--version"], [samtools_path]):
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=10,
            )
            version_text = "\n".join(
                part for part in (result.stdout.strip(), result.stderr.strip()) if part
            )
            match = re.search(
                r"(?:samtools\s+|Version:\s*)(\d+)\.(\d+)",
                version_text,
                flags=re.IGNORECASE,
            )
            if match:
                break
        log.debug(f"samtools version string: {version_text.splitlines()[:2]}")
        if not match:
            raise DependencyError(
                "Could not determine the samtools version; >= 1.12 is required."
            )
        major, minor = int(match.group(1)), int(match.group(2))
        if major < 1 or (major == 1 and minor < 12):
            raise DependencyError(
                f"samtools {major}.{minor} found, but >= 1.12 is required "
                "because mate rescue uses 'samtools view -N'."
            )
        log.debug(f"samtools {major}.{minor}: OK (>= 1.12).")
    except DependencyError:
        raise
    except Exception as e:
        raise DependencyError(
            f"samtools version check failed ({e}); >= 1.12 is required."
        ) from e


# ==============================================================================
# Command execution
# ==============================================================================


def run_command(cmd_list: list, log_file_handle, description: str = "") -> None:
    """Execute a shell command, raising PipelineStepError on failure.

    Does NOT call sys.exit, so the batch loop continues for other samples.

    Args:
        cmd_list:        Command and arguments as a list.
        log_file_handle: Open file handle to receive stdout and stderr.
        description:     Human-readable label for log messages.
    """
    label = f" ({description})" if description else ""
    rendered = shlex.join(str(part) for part in cmd_list)
    log.debug(f"CMD{label}: {rendered}")
    try:
        subprocess.run(
            cmd_list,
            stdout=log_file_handle, stderr=log_file_handle,
            text=True, check=True,
        )
        log.debug(f"CMD{label}: completed (rc=0)")
    except FileNotFoundError:
        raise PipelineStepError(f"Command not found: '{cmd_list[0]}'.")
    except subprocess.CalledProcessError as e:
        raise PipelineStepError(
            f"Command failed (rc={e.returncode}){label}: {rendered}"
        )
