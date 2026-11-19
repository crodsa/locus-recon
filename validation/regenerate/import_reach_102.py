#!/usr/bin/env python3
"""Which validation runners can reach ``reconcile_copy_number``.

Static analysis of the release tree, without running anything.  For every
runner the script follows (a) the imports inside ``locus_recon`` and
``validation`` and (b) the commands of the package that the reachable
validation code names (console scripts of pyproject.toml and ``-m`` modules),
and lists the reachable modules that refer to ``reconcile_copy_number``.  The
function is defined in ``locus_recon/copy_number.py``; release 1.0.2 changed
lines inside that function only (and the version string).

A runner that reaches no module referring to the function cannot be affected
by the change.  The analysis reads names in the source: a call built from a
string at run time would not be seen.

Usage: python import_reach_102.py <release tree> <output.tsv>
"""
import ast
import csv
import re
import sys
import tomllib
from pathlib import Path

FUNCTION = "reconcile_copy_number"
RUNNERS = [
    ("simulated benchmark", ["validation/run_mock_validation.py",
                             "validation/evaluate_mock_results.py"], "yes"),
    ("completeness benchmark", ["validation/run_completeness_benchmark.py"], "no"),
    ("tier calibration", ["validation/run_tier_calibration.py"], "no"),
    ("checks of single components", ["validation/run_frame_and_graph_validation.py"], "no"),
    ("synthetic copy series", ["validation/workflows/run_constructed_copy_series.py"], "no"),
    ("aphA1", ["validation/workflows/run_aphA1_validation.py"], "no"),
    ("23S rRNA gene and gyrB", ["validation/workflows/run_low_copy_23S.py"], "yes"),
    ("LIBA-6656 depth", ["validation/workflows/run_liba6656_depth.py"], "no"),
    ("LIBA-6656 overlap control", ["validation/depth-copy-number/overlap_null/*.py"], "no"),
    ("LIBA-6656 standard workflow", ["locus_recon/cli.py", "locus_recon/graph_cli.py",
                                     "validation/liba6656-gfa-rerun/verify_release_rerun.py"],
     "yes"),
    ("LIBA-6656 consensus without a graph", ["locus_recon/copy_number_cli.py"], "yes"),
    ("comparison with ARIBA and KMA", ["validation/run_comparator_benchmark.py"], "no"),
    ("block size of the bootstrap", ["validation/run_block_size_sensitivity.py"], "no"),
    ("external evaluation (other track)", ["validation/run_external_evaluation.py"],
     "not in this task"),
    ("targeted tests (other track)", ["validation/run_targeted_tests.py"], "not in this task"),
]


def module_name(root, path):
    parts = list(path.relative_to(root).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def main():
    root, output = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
    files = {}
    for package in ("locus_recon", "validation"):
        for path in sorted((root / package).rglob("*.py")):
            files[module_name(root, path)] = path
    by_stem = {}
    for name in files:
        by_stem.setdefault(name.split(".")[-1], []).append(name)
    project = tomllib.loads((root / "pyproject.toml").read_text())
    scripts = {command: target.split(":")[0]
               for command, target in project.get("project", {}).get("scripts", {}).items()}

    imports, mentions, commands, named_scripts = {}, {}, {}, {}
    for name, path in files.items():
        source = path.read_text()
        tree = ast.parse(source)
        found = set()
        package = name if path.name == "__init__.py" else name.rpartition(".")[0]
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    found.add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                base = node.module or ""
                if node.level:
                    parent = package.split(".")
                    parent = parent[:len(parent) - (node.level - 1)]
                    base = ".".join(parent + ([base] if base else []))
                found.add(base)
                for alias in node.names:
                    found.add(f"{base}.{alias.name}")
        resolved = set()
        for target in found:
            if target in files:
                resolved.add(target)
            # scripts run from validation/ import their neighbours by bare name
            elif "." not in target and len(by_stem.get(target, [])) == 1 and \
                    by_stem[target][0].startswith("validation"):
                resolved.add(by_stem[target][0])
        # a validation script that names another validation script (to start it,
        # or only in a comment) is taken to reach it: this errs on the side of
        # reporting a path that does not exist
        named = set()
        if name.startswith("validation"):
            for other, other_path in files.items():
                if other != name and other.startswith("validation") and \
                        other_path.name != "__init__.py" and \
                        re.search(rf"(?<!\w){re.escape(other_path.name)}\b", source):
                    named.add(other)
        named_scripts[name] = named - resolved
        imports[name] = resolved
        mentions[name] = bool(re.search(rf"\b{FUNCTION}\b", source)) and \
            name != "locus_recon.copy_number"
        started = set()
        for command, target in scripts.items():
            if re.search(rf"(?<![\w-]){re.escape(command)}(?![\w-])", source):
                started.add(target)
        for target in re.findall(r"[\"']-m[\"'],\s*[\"'](locus_recon\.[\w.]+)[\"']", source):
            started.add(target)
        for target in re.findall(r"-m\s+(locus_recon\.[\w.]+)", source):
            started.add(target)
        commands[name] = {target for target in started if target in files}

    rows = []
    for analysis, patterns, rerun in RUNNERS:
        starts = []
        for pattern in patterns:
            starts += [module_name(root, path) for path in sorted(root.glob(pattern))]
        missing = [pattern for pattern in patterns if not list(root.glob(pattern))]
        def closure(follow_named):
            reach, queue, started = set(), list(starts), set()
            while queue:
                name = queue.pop()
                if name in reach or name not in files:
                    continue
                reach.add(name)
                queue.extend(imports[name])
                if follow_named:
                    queue.extend(named_scripts[name])
                # only validation code and command-line modules start commands
                if name.startswith("validation") or name in starts:
                    started |= commands[name]
                    queue.extend(commands[name])
            return reach, started

        reach, started = closure(False)
        referring = sorted(name for name in reach if mentions[name])
        wide, _ = closure(True)
        only_named = sorted(name for name in wide if mentions[name] and name not in referring)
        via = sorted(f"{name} names {other.rsplit('.', 1)[-1]}.py"
                     for name in wide for other in named_scripts[name] if other in wide)
        rows.append([
            analysis, "; ".join(patterns) + (f" (not found: {'; '.join(missing)})" if missing else ""),
            "; ".join(sorted(started)) or "none",
            "yes" if "locus_recon.copy_number" in reach else "no",
            "; ".join(referring) or "none",
            "yes" if referring else "no",
            "; ".join(only_named) or "none",
            "; ".join(via) if only_named else "",
            rerun])
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["analysis", "runner", "package_commands_started",
                         "module_copy_number_reachable_by_imports_and_commands",
                         f"modules_reached_by_imports_and_commands_that_refer_to_{FUNCTION}",
                         f"can_reach_{FUNCTION}_by_imports_and_commands",
                         "further_referring_modules_if_scripts_named_in_text_are_followed",
                         "scripts_named_in_text",
                         "rerun_with_1.0.2"])
        writer.writerows(rows)
    for row in rows:
        print(" | ".join(row))
    print("console scripts:", scripts)


if __name__ == "__main__":
    main()
