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

    v1.0 used check=True with no fallback. The released source tarballs ship
    the full `validation/` tree but no `.git`, so every workflow raised here
    when run from the distributed artifact -- after the results had been
    written but before the provenance record was, leaving a reviewer with a
    non-zero exit and no provenance.json for the run they had just completed.

    Absence of a commit is now recorded rather than raised, the convention
    `tool_version` already uses for a missing tool. A release that wants the
    commit preserved outside a checkout should write it to a COMMIT file at
    the repository root; that file is read here in preference to failing.
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


def build_provenance(
    *,
    workflow: str,
    repo_root: str | Path,
    inputs: Iterable[str | Path],
    parameters: dict,
    tools: dict[str, Sequence[str | os.PathLike[str]]],
) -> dict:
    input_paths = require_files(inputs)
    return {
        "schema_version": "1.0",
        "workflow": workflow,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(repo_root),
        "parameters": parameters,
        "inputs": [
            {
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in input_paths
        ],
        "tools": {name: command_version(command) for name, command in sorted(tools.items())},
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
