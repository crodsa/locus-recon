from types import SimpleNamespace

import pytest

from locus_recon.pipeline import build_spades_command
from locus_recon.utils import (
    DependencyError,
    _check_samtools_version,
)


def test_legacy_samtools_fallback_is_rejected(monkeypatch):
    responses = iter([
        SimpleNamespace(stdout="", stderr="unrecognized command '--version'"),
        SimpleNamespace(
            stdout="",
            stderr="Program: samtools\nVersion: 0.1.19\nUsage: samtools <command>",
        ),
    ])
    monkeypatch.setattr("locus_recon.utils.subprocess.run", lambda *a, **k: next(responses))
    with pytest.raises(DependencyError, match=">= 1.12"):
        _check_samtools_version("samtools")


def test_modern_samtools_version_passes(monkeypatch):
    response = SimpleNamespace(stdout="samtools 1.20\nUsing htslib 1.20", stderr="")
    monkeypatch.setattr("locus_recon.utils.subprocess.run", lambda *a, **k: response)
    _check_samtools_version("samtools")


def test_spades_command_is_version_independent_and_never_careful():
    command = build_spades_command(
        "spades.py", "out", 4, 16, r1="R1.fq.gz", r2="R2.fq.gz", singletons="S.fq.gz",
    )
    assert command == [
        "spades.py", "-s", "S.fq.gz", "-1", "R1.fq.gz", "-2", "R2.fq.gz",
        "-o", "out", "-t", "4", "-m", "16", "--phred-offset", "33",
    ]
    retry = build_spades_command(
        "spades.py", "out", 4, 16, r1="R1.fq.gz", r2="R2.fq.gz", single_cell=True,
    )
    assert retry[:2] == ["spades.py", "--sc"]
    assert "--careful" not in command + retry


def test_spades_command_with_singletons_only():
    command = build_spades_command("spades.py", "out", 2, 8, singletons="S.fq.gz")
    assert "-1" not in command and "-2" not in command
    assert command[1:3] == ["-s", "S.fq.gz"]
