"""Run-manifest helpers for reproducible research outputs."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from argparse import Namespace
from typing import Iterable, Tuple

from . import PROGRAM, VERSION


def sha256_file(path: str, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 checksum of a file without loading it into memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tool_version(name: str, path: str) -> str:
    candidates = {
        "blastn": [[path, "-version"]],
        "makeblastdb": [[path, "-version"]],
        "samtools": [[path, "--version"], [path]],
        "spades.py": [[path, "--version"]],
        "bwa-mem2": [[path, "version"], [path]],
        "bwa": [[path]],
        "tqdm": [[path, "--version"]],
    }.get(name, [[path, "--version"]])

    for command in candidates:
        try:
            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
            output = (completed.stdout or completed.stderr).strip()
            if output and "unrecognized command" not in output.lower():
                return output.splitlines()[0][:500]
        except (OSError, subprocess.SubprocessError):
            continue
    return "unknown"


def _file_stat(path: str) -> dict:
    try:
        stat = os.stat(path)
        return {
            "path": path,
            "size_bytes": stat.st_size,
            "mtime_epoch": round(stat.st_mtime, 3),
        }
    except OSError:
        return {"path": path, "missing": True}


def write_run_manifest(
    manifest_path: str,
    args: Namespace,
    samples: Iterable[Tuple[str, str, str, str]],
    tools: dict,
    bait_profile: dict,
    run_date: str,
) -> None:
    """Write parameters, inputs, tool versions, and bait checksum as JSON."""
    sample_records = []
    for sample_id, assembly, r1, r2 in samples:
        sample_records.append({
            "sample_id": sample_id,
            "assembly": _file_stat(assembly),
            "r1": _file_stat(r1),
            "r2": _file_stat(r2),
        })

    serializable_profile = {
        key: value for key, value in bait_profile.items()
        if key not in {"lengths"}
    }
    manifest = {
        "schema_version": "1.0",
        "program": PROGRAM,
        "version": VERSION,
        "started_at": run_date,
        "command": sys.argv,
        "platform": {
            "python": platform.python_version(),
            "system": platform.platform(),
        },
        "parameters": vars(args),
        "bait": {
            **_file_stat(args.bait),
            "sha256": sha256_file(args.bait),
            "profile": serializable_profile,
        },
        "samplesheet": {
            **_file_stat(args.samplesheet),
            "sha256": sha256_file(args.samplesheet),
        },
        "samples": sample_records,
        "tools": {
            name: {"path": path, "version": _tool_version(name, path)}
            for name, path in sorted(tools.items())
            if path
        },
    }
    with open(manifest_path, "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
