from types import SimpleNamespace

import pytest

from locus_recon.utils import (
    DependencyError,
    _check_samtools_version,
    _check_spades_version,
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


def test_spades_4_disables_removed_careful_option(monkeypatch):
    response = SimpleNamespace(stdout="SPAdes genome assembler v4.2.0", stderr="")
    monkeypatch.setattr("locus_recon.utils.subprocess.run", lambda *a, **k: response)
    assert not _check_spades_version("spades.py")
