#!/usr/bin/env python3
"""DEPRECATED — superseded by the Java `eu.solven.matmul.docs.verify.CompareReferenceCatalogs`,
which compares against the UNION of FMM-Lille AND Perminov (this script is FMM-only),
reuses the manifest's field-correct vs-both flags, and emits the cited-vs-derived
classification. Kept only as a quick standalone cross-check. Prefer:
    mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.CompareReferenceCatalogs

Cross-reference our schemes/ catalog against fmm-lille's best-rank table
and report discrepancies.

For each format (n, m, p) sorted lexicographically:
- "ours"      = lowest rank among our verified schemes for that format
- "fmm"       = fmm-lille's best-known rank (from `references/catalogs/fmm-lille-catalog.json`)
- "delta"     = ours − fmm (positive means we're sub-optimal)

Writes `references/fmm-lille-discrepancies.md` and prints a summary.

Usage:
    python3 tools/compare_fmm_lille.py
"""
from __future__ import annotations

import json
import os
import re
import collections
from collections import defaultdict

CATALOG_PATH = "references/catalogs/fmm-lille-catalog.json"
SCHEMES_DIR = "src/main/resources/schemes"
OUT = "references/fmm-lille-discrepancies.md"

# Source may contain '_'; the shape is preceded by '-' (post-#173 rename) or
# legacy '_'. Rank token is _m{N} (or legacy _r{N}); a trailing _a{N}/field
# suffix may follow.
NAME_RE = re.compile(r"(?P<source>.*)[-_](?P<n>\d+)x(?P<m>\d+)x(?P<p>\d+)_(?:r|m)(?P<rank>\d+)[^/]*\.json")


def canonical(n: int, m: int, p: int) -> tuple[int, int, int]:
    return tuple(sorted((n, m, p)))


def load_fmm() -> dict[tuple[int, int, int], dict]:
    if not os.path.exists(CATALOG_PATH):
        raise FileNotFoundError(CATALOG_PATH)
    data = json.load(open(CATALOG_PATH))
    out = {}
    for row in data["entries"]:
        k = canonical(*row["format"])
        if k not in out or row["rank"] < out[k]["rank"]:
            out[k] = row
    return out


# fmm-lille's best-rank table is NON-commutative over a characteristic-0 field
# (Q). A scheme is comparable only if it is (a) non-commutative AND (b) valid
# over some characteristic-0 field. Excluding F₂/F₃-only schemes is load-bearing:
# e.g. AlphaTensor ⟨4,4,4⟩=47 is F₂-only, and its Kronecker power ⟨16³⟩=47²=2209
# would otherwise spuriously "beat" FMM's Q rank 2304 (field discipline).
CHAR0_FIELDS = {"Z", "Q", "R", "C"}


def comparable_to_fmm(path: str) -> bool:
    try:
        with open(path) as fh:
            j = json.load(fh)
    except (OSError, ValueError):
        return True  # can't read → don't silently drop; treat as comparable
    if j.get("commutative", False) is True:
        return False
    fields = j.get("fields")
    if fields is None:
        return True  # untagged → assume char-0 (legacy files default to Q/R/Z)
    return any(f in CHAR0_FIELDS for f in fields)


def is_atom(path: str) -> bool:
    """A scheme is an *atom* (primitive) iff its lineage root is absent or a
    single Atom/Leaf node — an explicit import or a formula-constructor ref
    (e.g. DIS09Lemma4(n=…)) — not composed by us from other catalog entries.
    Mirrors GenerateCatalogManifest.isAtomLineage. Lets the comparison segment
    'what we imported/derived-by-formula' from 'what we composed'."""
    try:
        with open(path) as fh:
            lin = json.load(fh).get("lineage")
    except (OSError, ValueError):
        return True
    if not isinstance(lin, dict):
        return True
    return lin.get("op", "") in ("", "Atom", "Leaf")


def load_ours() -> tuple[dict[tuple[int, int, int], list[dict]], int]:
    out: dict[tuple[int, int, int], list[dict]] = defaultdict(list)
    skipped = 0
    for root, _, files in os.walk(SCHEMES_DIR):
        for fn in files:
            if not fn.endswith(".json"):
                continue
            m = NAME_RE.match(fn)
            if not m:
                continue
            full = os.path.join(root, fn)
            if not comparable_to_fmm(full):
                skipped += 1
                continue
            n, mm, p = int(m["n"]), int(m["m"]), int(m["p"])
            k = canonical(n, mm, p)
            out[k].append({
                "source": m["source"],
                "rank": int(m["rank"]),
                "format": (n, mm, p),
                "file": os.path.relpath(full, SCHEMES_DIR),
                "atom": is_atom(full),
            })
    return out, skipped


def main() -> None:
    fmm = load_fmm()
    ours, skipped = load_ours()

    rows = []
    # Segmentation: (band, atom|composed) → Counter(tie/better/worse). Answers
    # "where do our FMM losses live" — they concentrate in *composed* entries
    # above MATERIALISE_MAX_DIM=16; the ≤16 band is atom-dominated and solved.
    seg = defaultdict(lambda: collections.Counter())
    for k, fmm_row in fmm.items():
        our_entries = ours.get(k, [])
        best_entry = min(our_entries, key=lambda e: e["rank"], default=None)
        our_best = best_entry["rank"] if best_entry else None
        delta = (our_best - fmm_row["rank"]) if our_best is not None else None
        if best_entry is not None:
            band = "<=16" if max(k) <= 16 else "17-32" if max(k) <= 32 else ">32"
            kind = "atom" if best_entry["atom"] else "composed"
            cmp = "tie" if delta == 0 else "better" if delta < 0 else "worse"
            seg[(band, kind)][cmp] += 1
        rows.append({
            "format": k,
            "max_dim": max(k),
            "fmm_rank": fmm_row["rank"],
            "fmm_refs": fmm_row.get("references", []),
            "our_rank": our_best,
            "our_count": len(our_entries),
            "delta": delta,
        })

    sub_optimal = [r for r in rows if r["delta"] is not None and r["delta"] > 0]
    matched = [r for r in rows if r["delta"] == 0]
    better = [r for r in rows if r["delta"] is not None and r["delta"] < 0]
    missing = [r for r in rows if r["our_rank"] is None]

    sub_optimal.sort(key=lambda r: (-r["delta"], r["format"]))
    better.sort(key=lambda r: (r["delta"], r["format"]))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write("# Catalog discrepancies vs fmm-lille\n\n")
        f.write(f"Generated from `{CATALOG_PATH}` (main fmm-lille catalog page).\n")
        f.write(f"Compares {sum(len(v) for v in ours.values())} of our scheme files "
                f"against fmm-lille's best-rank table.\n\n")
        f.write(f"- **Sub-optimal** (we have a higher rank than fmm-lille): {len(sub_optimal)}\n")
        f.write(f"- **Matched** (same rank): {len(matched)}\n")
        f.write(f"- **Better** (we have a lower rank than fmm-lille): {len(better)}\n")
        f.write(f"- **Missing** (fmm-lille has it; we don't): {len(missing)}\n\n")

        f.write("## Segmentation by atom vs composed (#198)\n\n")
        f.write("`atom` = our best entry is a primitive (explicit import or a "
                "formula-constructor ref); `composed` = we built it by Kron/concat/"
                "recombination/etc. Losses concentrate in *composed* entries above "
                "`MATERIALISE_MAX_DIM=16` — the ≤16 band is atom-dominated and solved.\n\n")
        f.write("| band | our entry | tie | better | worse |\n")
        f.write("|---|---|---:|---:|---:|\n")
        for band in ("<=16", "17-32", ">32"):
            for kind in ("atom", "composed"):
                c = seg[(band, kind)]
                if not c:
                    continue
                f.write(f"| {band} | {kind} | {c['tie']} | {c['better']} | {c['worse']} |\n")
        f.write("\n")

        f.write("## Sub-optimal (top 50)\n\n")
        f.write("Formats where fmm-lille knows of a strictly lower-rank algorithm than what we have.\n\n")
        f.write("| format | our rank | fmm rank | Δ | our count | fmm refs | fmm detail |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in sub_optimal[:50]:
            n, m, p = r["format"]
            refs = ", ".join(r["fmm_refs"][:3]) or "—"
            url = f"https://fmm.univ-lille.fr/{n}x{m}x{p}.html"
            f.write(f"| ⟨{n},{m},{p}⟩ | {r['our_rank']} | {r['fmm_rank']} | +{r['delta']} | {r['our_count']} | {refs} | [link]({url}) |\n")
        f.write("\n")

        if better:
            f.write("## Cases where we have a lower rank than fmm-lille\n\n")
            f.write("Likely indicates fmm-lille is missing a recent result — worth flagging upstream.\n\n")
            f.write("| format | our rank | fmm rank | Δ | our count |\n")
            f.write("|---|---|---|---|---|\n")
            for r in better[:30]:
                n, m, p = r["format"]
                f.write(f"| ⟨{n},{m},{p}⟩ | {r['our_rank']} | {r['fmm_rank']} | {r['delta']} | {r['our_count']} |\n")
            f.write("\n")

        if missing:
            missing.sort(key=lambda r: (r["max_dim"], r["format"]))
            f.write(f"## Missing from our catalog (top 30 / {len(missing)})\n\n")
            f.write("Formats in fmm-lille that we don't have any scheme for. "
                    "Ordered by max-dim ascending.\n\n")
            f.write("| format | fmm rank | fmm refs |\n")
            f.write("|---|---|---|\n")
            for r in missing[:30]:
                n, m, p = r["format"]
                refs = ", ".join(r["fmm_refs"][:3]) or "—"
                f.write(f"| ⟨{n},{m},{p}⟩ | {r['fmm_rank']} | {refs} |\n")

    print(f"sub-optimal={len(sub_optimal)} matched={len(matched)} "
          f"better={len(better)} missing={len(missing)} "
          f"(skipped {skipped} non-comparable files (commutative or non-char0) — fmm-lille is non-commutative)")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
