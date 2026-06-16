#!/usr/bin/env python3
"""
Flag synthesized / unverified entries in references/fmm-lille-catalog.json.

The digest mixes genuinely-scraped FMM ranks with auto-computed bounds. Two
tells, confirmed 2026-06-04:
  * the per-page detail scrape is EMPTY (raw_len==0) -> no FMM scheme there, so
    a rank present is fabricated (e.g. broken page <2,10,15>);
  * the rank equals the Hopcroft-Kerr closed form  ceil((3xy+max(x,y))/2)  for a
    <2,x,y> shape, with NO backing scheme realising it -> a formula bound, not a
    construction (e.g. <2,10,16>=248 while the real FMM scheme is 249).

We DON'T delete (the bound has provenance value) -- we annotate each entry with
  "verified": bool          # false = synthesized / not backed by a real scheme
  "rank_source": str        # why
so downstream (GenerateFmmGapReport) can refuse to treat unverified bounds as
real FMM results. naive_rank is a bogus copy of rank across all entries -> drop.
"""
import glob
import json
import os
import re

CAT = "references/fmm-lille-catalog.json"
DET = "references/fmm-lille-detail.json"


def hk_formula(fmt):
    """HK closed form for a <2,x,y> shape (the 2 on any axis), else None."""
    f = list(fmt)
    if 2 not in f:
        return None
    f.remove(2)
    if len(f) != 2:
        return None
    x, y = f
    return (3 * x * y + max(x, y) + 1) // 2


def main():
    cat = json.load(open(CAT))
    det = json.load(open(DET))
    det_items = det["entries"] if isinstance(det, dict) and "entries" in det else det
    rawlen = {x["url"]: x.get("raw_len", 0) for x in det_items}

    # best on-disk rank per shape across ALL catalog scheme files (encodes _m{rank}_).
    ondisk = {}
    for fpath in glob.glob("src/main/resources/schemes/**/*.json", recursive=True):
        m = re.search(r"-(\d+)x(\d+)x(\d+)_m(\d+)", os.path.basename(fpath))
        if not m:
            continue
        n, mm, p, r = map(int, m.groups())
        # consider all 6 axis permutations as the same shape (matmul symmetry)
        key = tuple(sorted((n, mm, p)))
        ondisk[key] = min(r, ondisk.get(key, 10 ** 9))

    flagged = 0
    by_reason = {}
    for e in cat["entries"]:
        fmt = e["format"]
        r = e["rank"]
        e.pop("naive_rank", None)  # bogus (== rank for all 5426 entries)
        empty = rawlen.get(e.get("details_url"), 0) == 0
        hk = hk_formula(fmt)
        is_formula = (hk is not None and r == hk)
        best = ondisk.get(tuple(sorted(fmt)))
        backed = best is not None and best <= r  # some real scheme is this good or better
        unverified = (empty or is_formula) and not backed
        if unverified:
            e["verified"] = False
            reason = "empty-page" if empty else "hk-formula"
            if empty and is_formula:
                reason = "empty-page+hk-formula"
            e["rank_source"] = "synthesized:" + reason
            flagged += 1
            by_reason[reason] = by_reason.get(reason, 0) + 1
        else:
            e["verified"] = True

    json.dump(cat, open(CAT, "w"), indent=1)
    total = len(cat["entries"])
    print("entries: %d   flagged unverified: %d   (verified: %d)" % (total, flagged, total - flagged))
    print("by reason:", by_reason)


if __name__ == "__main__":
    main()
