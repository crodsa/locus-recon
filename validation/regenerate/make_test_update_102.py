#!/usr/bin/env python3
"""Suggested update of one test for the 23S records of release 1.0.2.

``test_committed_low_copy_panel_meets_prespecified_two_copy_boundary`` reads
``validation/low-copy-23S/results.tsv`` and ``evaluation.json``.  With release
1.0.2 the four graph counts of two are returned for review, so the runner no
longer counts them as exact supported calls (``exact_two_copy_consensus_calls``
is 0, not 4).  The update changes that expected value and adds assertions on
the new status and bounds of the four loci; no assertion is removed.

Run from the root of the release tree; it edits the test file in place, and
``git diff -- tests`` is then the patch.
"""
import re
from pathlib import Path

NAME = "test_committed_low_copy_panel_meets_prespecified_two_copy_boundary"
CHANGES = [
    ('"""The depth criteria hold; one 23S graph count is only a lower bound.',
     '"""The depth criteria hold; no 23S graph count is returned as supported.'),
    ("so the count is reported as a lower bound and the exact-call criterion is "
     "met in four of five loci. No reported integer or bound contradicts truth.",
     "so the count is reported as a lower bound. At the other four loci the graph\n"
     "    holds two contexts and the depth-ratio interval excludes two, so the count\n"
     "    is returned for review with a lower bound of two and no upper bound.\n"
     "    No reported integer or bound contradicts truth."),
    ('assert all(row["copy_number_method"] == "GRAPH_COUNT_OVER_DISCORDANT_DEPTH" '
     "for row in exact)",
     'assert all(row["copy_number_method"] == "GRAPH_COUNT_OVER_DISCORDANT_DEPTH"\n'
     "               for row in exact)\n"
     '    assert all(row["copy_number_status"] == "REVIEW" for row in exact)\n'
     '    assert all((row["copy_number_lower_bound"], row["copy_number_upper_bound"]) == ("2", "")\n'
     "               for row in exact)\n"
     '    assert all(int(row["copy_number_lower_bound"]) <= int(row["truth_copies"])\n'
     "               for row in exact)"),
    ('assert evaluation["exact_two_copy_consensus_calls"] == 4',
     'assert evaluation["exact_two_copy_consensus_calls"] == 0'),
]


def main():
    path = Path("tests/test_validation_workflows.py")
    text = path.read_text()
    start = text.index(f"def {NAME}(")
    end = text.find("\n\n\ndef ", start)
    end = len(text) if end == -1 else end
    body = text[start:end]
    for old, new in CHANGES:
        pattern = r"\s+".join(re.escape(word) for word in old.split())
        found = re.findall(pattern, body)
        if len(found) != 1:
            raise SystemExit(f"passage found {len(found)} times in {NAME}: {old[:60]}")
        body = re.sub(pattern, lambda _: new, body)
    path.write_text(text[:start] + body + text[end:])
    print(f"{NAME}: {len(CHANGES)} passages changed")


if __name__ == "__main__":
    main()
