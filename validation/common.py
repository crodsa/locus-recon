"""Shared fail-fast and provenance helpers for validation workflows."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_files(paths: Iterable[str | Path]) -> list[Path]:
    resolved = [Path(path).resolve() for path in paths]
    missing = [str(path) for path in resolved if not path.is_file()]
    if missing:
        raise FileNotFoundError("required input file(s) missing:\n  " + "\n  ".join(missing))
    return resolved


def run_checked(
    command: Sequence[str | os.PathLike[str]],
    *,
    cwd: str | Path | None = None,
    stdout=None,
    stderr=None,
) -> subprocess.CompletedProcess:
    cmd = [str(part) for part in command]
    completed = subprocess.run(cmd, cwd=cwd, stdout=stdout, stderr=stderr, text=True)
    if completed.returncode != 0:
        detail = completed.stderr if isinstance(completed.stderr, str) else ""
        raise RuntimeError(
            f"command failed with exit {completed.returncode}: {' '.join(cmd)}"
            + (f"\n{detail[-4000:]}" if detail else "")
        )
    return completed


def command_version(command: Sequence[str | os.PathLike[str]]) -> dict:
    cmd = [str(part) for part in command]
    try:
        completed = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"command": cmd, "version": "UNAVAILABLE", "error": str(exc)}
    text = (completed.stdout or completed.stderr).strip()
    return {
        "command": cmd,
        "version": text.splitlines()[0] if text else f"exit={completed.returncode}",
        "exit_code": completed.returncode,
    }


def git_commit(repo_root: str | Path) -> str:
    """Commit of the working tree, or a recorded/UNAVAILABLE sentinel.

    Source archives ship ``validation/`` without ``.git``, so the absence of a
    commit is recorded rather than raised: a workflow must still write its
    provenance after writing its results.  A COMMIT file at the repository
    root, when present, is read in preference to the sentinel.
    """
    recorded = Path(repo_root) / "COMMIT"
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo_root,
            capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        completed = None
        reason = str(exc)
    else:
        if completed.returncode == 0 and completed.stdout.strip():
            return completed.stdout.strip()
        reason = (completed.stderr or "").strip().splitlines()[0:1] or ["no HEAD"]
        reason = reason[0]
    if recorded.is_file():
        text = recorded.read_text().strip().split()
        if text:
            return f"{text[0]} (recorded; not a git checkout)"
    return f"UNAVAILABLE ({reason})"


def source_sha256(repo_root: str | Path) -> str:
    """SHA-256 over the package sources, independent of git.

    Every ``locus_recon/*.py`` file contributes its repository-relative path
    and its bytes, in sorted order, so the digest identifies the code that ran
    whether or not it was committed, and anyone holding a release can
    recompute it.
    """
    root = Path(repo_root)
    digest = hashlib.sha256()
    for path in sorted((root / "locus_recon").rglob("*.py")):
        digest.update(path.relative_to(root).as_posix().encode("utf-8") + b"\0")
        digest.update(path.read_bytes() + b"\0")
    return digest.hexdigest()


def software_record(repo_root: str | Path) -> dict:
    """Package version, source digest and git state of the code that ran."""
    root = Path(repo_root)
    version = "unknown"
    init = root / "locus_recon" / "__init__.py"
    if init.is_file():
        for line in init.read_text().splitlines():
            if line.startswith("VERSION ="):
                version = line.split("=", 1)[1].strip().strip("\"'")
                break
    record = {
        "name": "locus-recon",
        "version": version,
        "source_sha256": source_sha256(root) if (root / "locus_recon").is_dir() else None,
        "git_commit": git_commit(root),
        "git_tree_clean": None,
    }
    try:
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no", "--",
             "locus_recon", "validation"],
            cwd=root, capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return record
    if status.returncode == 0:
        record["git_tree_clean"] = not status.stdout.strip()
    return record


def portable_path(path: str | Path, repo_root: str | Path) -> str:
    """Repository-relative path, or the file name for inputs outside it.

    Absolute paths describe the machine a run happened on, not the evidence;
    the checksum beside each entry is what identifies an input.
    """
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(Path(repo_root).resolve()).as_posix()
    except ValueError:
        return resolved.name


def build_provenance(
    *,
    workflow: str,
    repo_root: str | Path,
    inputs: Iterable[str | Path],
    parameters: dict,
    tools: dict[str, Sequence[str | os.PathLike[str]]],
) -> dict:
    input_paths = require_files(inputs)
    software = software_record(repo_root)
    return {
        "schema_version": "1.0",
        "workflow": workflow,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": software["git_commit"],
        "software": software,
        "parameters": parameters,
        "inputs": [
            {
                "path": portable_path(path, repo_root),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in input_paths
        ],
        "tools": {
            name: portable_command_version(command, repo_root)
            for name, command in sorted(tools.items())
        },
    }


def portable_command_version(
    command: Sequence[str | os.PathLike[str]], repo_root: str | Path,
) -> dict:
    """``command_version`` with an absolute executable reduced to its name."""
    record = command_version(command)
    record["command"] = [
        Path(part).name if os.path.isabs(part) else part for part in record["command"]
    ]
    return record


def compact_summary(provenance: dict, **extra) -> dict:
    """Short run record for a deposit: what ran, with which code and tools.

    The full provenance, with every input checksum, stays with the run's own
    output directory; the deposit keeps this summary beside the results.
    """
    software = provenance["software"]
    return {
        "schema_version": "1.0",
        "workflow": provenance["workflow"],
        "created_utc": provenance["created_utc"],
        "software": {key: software[key] for key in ("name", "version", "source_sha256")},
        "parameters": provenance["parameters"],
        "tools": {name: record["version"] for name, record in provenance["tools"].items()},
        **extra,
    }


def write_json(path: str | Path, payload: object) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=destination.parent, delete=False
    ) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(destination)
