#!/usr/bin/env python3
"""For each format in fmm-lille's catalog, compute our best achievable rank
via Kronecker composition of our existing catalog schemes, then categorize:

- `derive_wins`   — our derived/composed best beats fmm-lille's listed rank
- `derive_match`  — same rank as fmm-lille
- `fmm_wins`      — fmm-lille has a strictly lower rank (these are the
                    legitimate import candidates)
- `no_coverage`   — neither we nor fmm has anything (won't appear unless
                    we ask for an out-of-fmm format)

Methodology — dynamic programming over canonical (sorted) formats:
    best[n,m,p] = min(direct catalog ranks for any ordering)
                ∪ {best[a,b,c] · best[n/a, m/b, p/c]
                   for every factorization with both sides known}

Iterates to fixed point.

Output: `references/best-composition-report.md` and a list of
"import-from-fmm" candidates that respects "fmm only if better".

Usage:
    python3 tools/find_best_compositions.py
"""
from __future__ import annotations

import json
import os
import re
from collections import defaultdict

SCHEMES_DIR = "src/main/resources/schemes"
FMM_CATALOG = "references/catalogs/fmm-lille-catalog.json"
OUT_REPORT = "references/best-composition-report.md"
OUT_IMPORT_LIST = "references/import-from-fmm.txt"

# Match new-form filenames: <source>_<n>x<m>x<p>_r<rank>_a<adds>[_tag].json
NAME_RE = re.compile(
    r"(?P<source>.*)_(?P<n>\d+)x(?P<m>\d+)x(?P<p>\d+)_(?:r|m)(?P<rank>\d+)(?:_a\d+)?[^/]*\.json"
)


def canonical(n: int, m: int, p: int) -> tuple[int, int, int]:
    return tuple(sorted((n, m, p)))


def load_direct_catalog() -> dict[tuple[int, int, int], int]:
    """Walk src/main/resources/schemes/, return canonical-format → best direct rank."""
    out: dict[tuple[int, int, int], int] = {}
    for root, _, files in os.walk(SCHEMES_DIR):
        for fn in files:
            if not fn.endswith(".json"):
                continue
            m = NAME_RE.match(fn)
            if not m:
                continue
            n, mm, p = int(m["n"]), int(m["m"]), int(m["p"])
            rank = int(m["rank"])
            k = canonical(n, mm, p)
            if k not in out or rank < out[k]:
                out[k] = rank
    # The trivial scalar case is always rank 1.
    out.setdefault((1, 1, 1), 1)
    return out


def load_fmm_catalog() -> dict[tuple[int, int, int], int]:
    if not os.path.exists(FMM_CATALOG):
        return {}
    data = json.load(open(FMM_CATALOG))
    out: dict[tuple[int, int, int], int] = {}
    for row in data["entries"]:
        k = canonical(*row["format"])
        if k not in out or row["rank"] < out[k]:
            out[k] = row["rank"]
    return out


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def fixed_point_compose(direct: dict, max_dim: int = 32):
    """DP over factorizations. Returns (best, recipe) per canonical format.

    recipe[k] = ("direct", source) or ("compose", left_format, right_format)
    """
    best = dict(direct)
    recipe: dict[tuple[int, int, int], tuple] = {k: ("direct",) for k in direct}

    iteration = 0
    while True:
        iteration += 1
        changed = 0
        # Enumerate every canonical target.
        for n in range(1, max_dim + 1):
            for m in range(n, max_dim + 1):
                for p in range(m, max_dim + 1):
                    target = (n, m, p)
                    # Try all ordered factorizations (a1·a2=n, b1·b2=m, c1·c2=p).
                    for a1 in divisors(n):
                        a2 = n // a1
                        for b1 in divisors(m):
                            b2 = m // b1
                            for c1 in divisors(p):
                                c2 = p // c1
                                f1 = canonical(a1, b1, c1)
                                f2 = canonical(a2, b2, c2)
                                if f1 == (1, 1, 1) and f2 == target:
                                    continue
                                if f2 == (1, 1, 1) and f1 == target:
                                    continue
                                if f1 not in best or f2 not in best:
                                    continue
                                cand = best[f1] * best[f2]
                                if cand < best.get(target, 10**18):
                                    best[target] = cand
                                    recipe[target] = ("compose", f1, f2)
                                    changed += 1
        if changed == 0:
            break
    return best, recipe


def write_report(direct, best, recipe, fmm) -> None:
    # Categorize per fmm format.
    derive_wins = []
    derive_match = []
    fmm_wins = []

    for fmt, fmm_rank in fmm.items():
        our = best.get(fmt)
        if our is None:
            fmm_wins.append((fmt, None, fmm_rank))
            continue
        if our < fmm_rank:
            derive_wins.append((fmt, our, fmm_rank))
        elif our == fmm_rank:
            derive_match.append((fmt, our, fmm_rank))
        else:
            fmm_wins.append((fmt, our, fmm_rank))

    derive_wins.sort(key=lambda r: (r[2] - r[1], r[0]), reverse=True)
    fmm_wins.sort(key=lambda r: ((r[1] or 10**18) - r[2], r[0]), reverse=True)

    with open(OUT_REPORT, "w") as f:
        f.write("# Best composed rank vs fmm-lille\n\n")
        f.write("For each fmm-lille format, we compute our best achievable rank via "
                "either a direct scheme on disk or a Kronecker composition of two "
                "smaller schemes (dynamic programming over all factorizations). "
                "Then categorize.\n\n")
        f.write(f"- **derive-wins** (our composed/direct beats fmm): {len(derive_wins)}\n")
        f.write(f"- **derive-match** (same rank): {len(derive_match)}\n")
        f.write(f"- **fmm-wins** (legitimate import candidates): {len(fmm_wins)}\n\n")

        f.write("## Derive-wins — formats where our catalog implies a lower rank than fmm-lille's main page lists\n\n")
        f.write("Likely fmm-lille hasn't propagated a recent composition.\n\n")
        f.write("| format | our rank | fmm rank | Δ | recipe |\n|---|---|---|---|---|\n")
        for fmt, our, fmm_rank in derive_wins[:40]:
            n, m, p = fmt
            r = recipe.get(fmt, ("?",))
            if r[0] == "compose":
                f1, f2 = r[1], r[2]
                rec = f"⟨{f1[0]},{f1[1]},{f1[2]}⟩ ⊗ ⟨{f2[0]},{f2[1]},{f2[2]}⟩"
            else:
                rec = "direct"
            f.write(f"| ⟨{n},{m},{p}⟩ | {our} | {fmm_rank} | -{fmm_rank - our} | {rec} |\n")

        f.write("\n## fmm-wins — legitimate import candidates (fmm has strictly lower rank than we can derive)\n\n")
        f.write("These are the formats worth importing from fmm-lille's `.mpl.bz2`. "
                "After running `tools/import_fmm_maple.py {n}x{m}x{p}` for each, "
                "the discrepancy report should empty out.\n\n")
        f.write("| format | our rank | fmm rank | Δ |\n|---|---|---|---|\n")
        for fmt, our, fmm_rank in fmm_wins[:80]:
            n, m, p = fmt
            our_txt = "—" if our is None else our
            delta = "—" if our is None else f"+{our - fmm_rank}"
            f.write(f"| ⟨{n},{m},{p}⟩ | {our_txt} | {fmm_rank} | {delta} |\n")

    # Plain-text import list for piping into the importer.
    with open(OUT_IMPORT_LIST, "w") as f:
        for fmt, _, _ in fmm_wins:
            f.write(f"{fmt[0]}x{fmt[1]}x{fmt[2]}\n")

    print(f"derive-wins={len(derive_wins)} derive-match={len(derive_match)} fmm-wins={len(fmm_wins)}")
    print(f"wrote {OUT_REPORT}")
    print(f"wrote {OUT_IMPORT_LIST}")


def main():
    direct = load_direct_catalog()
    print(f"direct catalog: {len(direct)} canonical formats")
    fmm = load_fmm_catalog()
    print(f"fmm catalog: {len(fmm)} canonical formats")

    best, recipe = fixed_point_compose(direct, max_dim=32)
    print(f"best after composition DP: {len(best)} canonical formats")

    write_report(direct, best, recipe, fmm)


if __name__ == "__main__":
    main()
