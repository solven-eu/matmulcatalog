#!/usr/bin/env python3
"""
Filename rename per task #173: source-to-shape uses `-`, source-internal
uses `_`, drop trailing-field tokens.

Strategy: locate the SHAPE token `_NxMxP_(m|r)N` in each filename. Split
into (source, shape_onward). In `source`: replace every `-` with `_`.
Stitch back with `-` between source and shape. In `shape_onward`: strip
any TRAILING `_<field>` token after the rank/additions (info now lives
in JSON `fields[]` per #174).

Trailing field tokens to strip: `_Q`, `_R`, `_C`, `_Z`, `_ZT`, `_F2`,
`_F3`, `_0.5xC`, `_0.5xZ`. Source-INTERNAL field tokens (before the
shape) stay as part of the source — e.g. `perminov-ZT_3x7x7` becomes
`perminov_ZT-3x7x7`.
"""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path("src/main/resources/schemes")
SHAPE_RE = re.compile(r"^(.+?)[_-](\d+x\d+x\d+_(?:m|r)\d+(?:_a\d+)?)(.*?)(\.json)$")

TRAILING_FIELD_RE = re.compile(
    r"_(?:F2|F3|ZT|Q|R|C|Z|0\.5xC|0\.5xZ|complex|commutative)$"
)


def rename_one(filename: str) -> str | None:
    """Return the new filename, or None if the file shouldn't be renamed."""
    m = SHAPE_RE.match(filename)
    if not m:
        return None
    source, shape, trailing, ext = m.group(1), m.group(2), m.group(3), m.group(4)

    # Source: all `-` → `_`
    source_new = source.replace("-", "_")

    # Trailing: strip field tokens iteratively (some files have e.g.
    # `_Q_commutative` — strip both).
    while True:
        m2 = TRAILING_FIELD_RE.search(trailing)
        if not m2:
            break
        trailing = trailing[: m2.start()]

    new_name = f"{source_new}-{shape}{trailing}{ext}"
    return new_name if new_name != filename else None


def main():
    apply = "--apply" in sys.argv

    files = sorted(ROOT.rglob("*.json"))
    renames = []
    collisions = defaultdict(list)

    for f in files:
        new_name = rename_one(f.name)
        if new_name is None:
            continue
        new_path = f.parent / new_name
        renames.append((f, new_path))
        collisions[new_path].append(f)

    print(f"{len(files)} files scanned; {len(renames)} renames computed")

    # Detect collisions
    actual_collisions = {dst: srcs for dst, srcs in collisions.items() if len(srcs) > 1}
    if actual_collisions:
        print(f"\n!! {len(actual_collisions)} collisions detected:")
        for dst, srcs in list(actual_collisions.items())[:10]:
            print(f"   {dst.name}  ← {[s.name for s in srcs]}")
        if not apply:
            print("\nNot applying. Fix collisions first.")
            sys.exit(1)

    # Show samples
    print("\nSample renames (first 10):")
    for old, new in renames[:10]:
        print(f"  {old.name}\n  → {new.name}")
    print(f"\n  …{len(renames) - 10} more")

    if apply:
        for old, new in renames:
            old.rename(new)
        print(f"\nApplied: {len(renames)} files renamed")


if __name__ == "__main__":
    main()
