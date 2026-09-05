"""Geometry-matched null for the additive dosage estimand.

The random-window control shows what a single-copy query returns when the
assembler represents it as one interval. It does not answer the question the
tcdB case raises: what does the same estimand return for a single-copy query
that the discovery step fragments into many mutually overlapping intervals,
which is the geometry that produces a dosage above one?

This script enriches for that geometry. It self-aligns the assembly once,
builds a per-base multiplicity track, and draws pseudo-queries centred on
repeat-rich positions, keeping those whose own discovery run returns at least
four intervals. Each retained query is single-copy in the sense that matters:
it is a stretch of the draft assembly, not a locus with independent evidence
for two copies. Any dosage above one that it returns is estimand behaviour,
not biology.
"""
import json, random, subprocess, sys, tempfile
from pathlib import Path
sys.path.insert(0, "lr_v11")
from locus_recon.depth_ratio import query_copy_profile, _query_segments

H = Path("/Users/crodsa/Library/CloudStorage/OneDrive-UniversidaddeCostaRica/Cesar/Papers/"
         "EnPreparacion/2026/Locus_Recon_LIBA6656_ST154/MG_technical_resource/"
         "Locus_Recon_MG_pre_submission_handoff_2026-08-31")
W = H / "02_DATA_AND_REPRODUCIBILITY/validation_working_data/LIBA6656/results"
ASM, BAM = W / "LIBA6656_ST154.assembly.fasta", W / "ERR467623_vs_LIBA6656_draft.sorted.bam"
QLEN, SEED, TARGET, MIN_REGIONS = 7104, 20260903, 400, 2
MAX_TRIED = 3000
LOCUS = {".11940_5_90.1", ".11940_5_90.29", ".11940_5_90.56", ".11940_5_90.59",
         ".11940_5_90.63", ".11940_5_90.76", ".11940_5_90.77", ".11940_5_90.8"}


def read_fasta(p):
    d, k = {}, None
    for ln in open(p):
        if ln.startswith(">"):
            k = ln[1:].split()[0]
            d[k] = []
        else:
            d[k].append(ln.strip())
    return {k: "".join(v) for k, v in d.items()}


def med(xs):
    s = sorted(xs)
    n = len(s)
    return 0.0 if not n else (s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2)


asm = read_fasta(ASM)
out = subprocess.run(["samtools", "depth", "-a", "-q", "20", "-Q", "20", str(BAM)],
                     capture_output=True, text=True, check=True)
dep = {}
for ln in out.stdout.splitlines():
    c, _, d = ln.split("\t")
    dep.setdefault(c, []).append(int(d))
bb = [d for c, v in dep.items() if c not in LOCUS and len(v) >= 5000 for d in v[100:-100]]
BB = med(bb)

db = tempfile.mkdtemp()
subprocess.run(["makeblastdb", "-in", str(ASM), "-dbtype", "nucl", "-out", f"{db}/asm"],
               check=True, capture_output=True)

# one self-alignment to locate repeat-rich positions
self_bl = subprocess.run(
    ["blastn", "-query", str(ASM), "-db", f"{db}/asm", "-evalue", "1e-10",
     "-num_threads", "4", "-outfmt", "6 qseqid sseqid pident length qstart qend"],
    capture_output=True, text=True, check=True)
mult = {c: [0] * len(s) for c, s in asm.items()}
for ln in self_bl.stdout.splitlines():
    f = ln.split("\t")
    if int(f[3]) < 300 or float(f[2]) < 90.0:
        continue
    qs, qe = int(f[4]) - 1, int(f[5])
    tr = mult.get(f[0])
    if tr is None:
        continue
    for i in range(qs, min(qe, len(tr))):
        tr[i] += 1

cands = []
for c, s in asm.items():
    if c in LOCUS or len(s) < QLEN + 400:
        continue
    tr = mult[c]
    for st in range(100, len(s) - QLEN - 100, 250):
        w = tr[st:st + QLEN]
        rich = sum(1 for x in w if x >= 2) / QLEN
        if s[st:st + QLEN].upper().count("N") == 0:
            cands.append((c, st, rich))
random.Random(SEED).shuffle(cands)
print("repeat-rich candidate windows:", len(cands), flush=True)

rows, tried = [], 0
q = Path(db) / "q.fa"
for c, st, rich in cands:
    if len(rows) >= TARGET or tried >= MAX_TRIED:
        break
    tried += 1
    q.write_text(">q\n%s\n" % asm[c][st:st + QLEN])
    bl = subprocess.run(["blastn", "-query", str(q), "-db", f"{db}/asm", "-evalue", "1e-10",
                         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send"],
                        capture_output=True, text=True, check=True)
    hits = []
    for ln in bl.stdout.splitlines():
        f = ln.split("\t")
        if int(f[3]) < 300 or float(f[2]) < 90.0:
            continue
        ss, se = int(f[6]), int(f[7])
        hits.append({"sseqid": f[1], "sstart": min(ss, se), "send": max(ss, se),
                     "qstart": int(f[4]), "qend": int(f[5])})
    merged, byc = [], {}
    for h in hits:
        byc.setdefault(h["sseqid"], []).append(h)
    for cc, hs in byc.items():
        hs.sort(key=lambda x: x["sstart"])
        cur = dict(hs[0])
        for h in hs[1:]:
            if h["sstart"] <= cur["send"] + 20:
                cur["send"] = max(cur["send"], h["send"])
                cur["qstart"] = min(cur["qstart"], h["qstart"])
                cur["qend"] = max(cur["qend"], h["qend"])
            else:
                merged.append(cur)
                cur = dict(h)
        merged.append(cur)
    if len(merged) < MIN_REGIONS:
        continue
    spans, ratios = [], []
    for h in merged:
        d = dep.get(h["sseqid"], [])[h["sstart"] - 1:h["send"]]
        if not d:
            continue
        spans.append((min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])))
        ratios.append(med(d) / BB)
    if len(spans) < MIN_REGIONS:
        continue
    segs = _query_segments(spans, QLEN)
    dose, lo, hi, covbp, mx = query_copy_profile(segs, ratios, QLEN)
    total = sum(e - s + 1 for s, e in spans)
    rows.append({"contig": c, "start": st, "repeat_fraction": round(rich, 3),
                 "n_regions": len(spans), "dosage": dose, "geom_lower": lo,
                 "geom_upper": hi, "covered_bp": covbp, "coverage_fraction": covbp / QLEN,
                 "overlap_bp": total - covbp,
                 "overlap_fraction": (total - covbp) / QLEN,
                 "max_multiplicity": mx})
    if len(rows) % 25 == 0:
        print("  kept %d / tried %d" % (len(rows), tried), flush=True)

json.dump({"backbone_median": BB, "query_len": QLEN, "seed": SEED,
           "min_regions": MIN_REGIONS, "windows_tried": tried,
           "n": len(rows), "rows": rows},
          open("out/null_dose_matched.json", "w"), indent=1)
print("written", len(rows))
