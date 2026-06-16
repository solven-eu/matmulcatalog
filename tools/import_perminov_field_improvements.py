"""
Bulk-import field-specific rank improvements from Perminov's status.json.

For each ⟨a,b,c⟩ shape, Perminov publishes separate ranks per field
(Z, Q, ZT). When the Q-arithmetic rank is strictly lower than what we
have on disk (typically because we only imported ZT), this script
chains to `import_fmm_maple.py <shape>` to fetch the upstream Maple
file — Perminov's `source` field for these Q entries always points to
`schemes/known/tensor/<shape>_tensor.mpl` which is FMM.

Strategy:
1. Load status.json (cached at /tmp/perminov_status.json).
2. For each shape × field, if Perminov rank < our best local rank:
   - if source endswith `_tensor.mpl` → run import_fmm_maple.py
   - else → log and skip (we don't have a fetcher for it yet)
3. Dry-run by default; pass --apply to actually invoke imports.

Usage:
  python tools/import_perminov_field_improvements.py            # dry-run
  python tools/import_perminov_field_improvements.py --apply    # do it
  python tools/import_perminov_field_improvements.py --max 5    # cap to 5
"""
import argparse
import json
import os
import re
import subprocess
import sys
import urllib.request


STATUS_URL = "https://raw.githubusercontent.com/dronperminov/FastMatrixMultiplication/master/schemes/status.json"
LOCAL_STATUS = "/tmp/perminov_status.json"
SCHEMES_ROOT = "src/main/resources/schemes"


def load_status(force_refresh: bool = False):
    """Cache Perminov status.json in /tmp. Re-fetch when force_refresh
    is set, or when the local copy is missing/empty.
    """
    if force_refresh or not os.path.exists(LOCAL_STATUS) or os.path.getsize(LOCAL_STATUS) < 1000:
        print(f"fetching {STATUS_URL}")
        urllib.request.urlretrieve(STATUS_URL, LOCAL_STATUS)
    with open(LOCAL_STATUS) as f:
        return json.load(f)


def load_local_ranks() -> dict[tuple[int, int, int], int]:
    """Lowest rank per (a,b,c) sorted-shape across all our scheme files."""
    pat = re.compile(r"_(\d+)x(\d+)x(\d+)_r(\d+)")
    out: dict[tuple[int, int, int], int] = {}
    for root, _, files in os.walk(SCHEMES_ROOT):
        for f in files:
            if not f.endswith(".json"):
                continue
            m = pat.search(f)
            if not m:
                continue
            n, a, b, r = (int(x) for x in m.groups())
            key = tuple(sorted([n, a, b]))
            if key not in out or r < out[key]:
                out[key] = r
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--apply", action="store_true",
                    help="actually invoke import scripts (default: dry-run)")
    ap.add_argument("--max", type=int, default=0,
                    help="cap total imports (0 = no cap)")
    ap.add_argument("--field", default="Q",
                    help="which Perminov field to pull (Q | Z | ZT, default Q)")
    ap.add_argument("--refresh-status", action="store_true",
                    help="re-fetch Perminov status.json (default: use cached if present)")
    args = ap.parse_args()

    status = load_status(force_refresh=args.refresh_status)
    local = load_local_ranks()

    candidates = []
    for shape_str, entry in status.items():
        m = re.match(r"^(\d+)x(\d+)x(\d+)$", shape_str)
        if not m:
            continue
        a, b, c = (int(x) for x in m.groups())
        shape = tuple(sorted([a, b, c]))
        ranks = entry.get("ranks", {})
        per_field = ranks.get(args.field)
        if per_field is None:
            continue
        local_best = local.get(shape)
        if local_best is not None and per_field >= local_best:
            continue
        # Look up source for that field's best scheme. If Perminov publishes
        # a rank in `ranks.<field>` without a corresponding `schemes.<field>[]`
        # entry, the rank typically reflects the FMM-published bound that
        # Perminov tracks but hasn't re-hosted. Fall back to FMM upstream.
        schemes = entry.get("schemes", {}).get(args.field, [])
        if schemes:
            src = schemes[0].get("source", "")
        else:
            src = f"schemes/known/tensor/{shape_str}_tensor.mpl"
        candidates.append((shape_str, shape, per_field, local_best, src))

    candidates.sort(key=lambda x: (x[1][2], x[1][1], x[1][0]))  # smallest first

    print(f"{'shape':<10} {'P.'+args.field:<6} {'ours':<6} source")
    print("-" * 78)
    for shape_str, shape, per_field, local_best, src in candidates:
        print(f"{shape_str:<10} {per_field:<6} {str(local_best):<6} {src}")
    print()
    print(f"Total improvements to fetch: {len(candidates)}")

    if not args.apply:
        print("\n(dry-run; re-run with --apply to fetch)")
        return 0

    done = 0
    for shape_str, shape, per_field, local_best, src in candidates:
        if args.max and done >= args.max:
            break
        if src.endswith("_tensor.mpl"):
            print(f"\n=== {shape_str}: importing from FMM Maple (Perminov points at upstream) ===")
            r = subprocess.run(
                ["/tmp/strassen_venv/bin/python", "tools/import_fmm_maple.py", shape_str],
                check=False, capture_output=True, text=True)
            print(r.stdout, end="")
            if r.returncode != 0:
                print(r.stderr, file=sys.stderr)
                print(f"  ERROR exit={r.returncode}")
                continue
            done += 1
        elif src.startswith("schemes/results/"):
            # Perminov's own scheme — fetch the JSON directly.
            print(f"\n=== {shape_str}: fetching Perminov's own JSON ({src}) ===")
            url = f"https://raw.githubusercontent.com/dronperminov/FastMatrixMultiplication/master/{src}"
            sec = max(shape)
            # Use _r{rank}_ naming (project convention) — keeps load_local_ranks()
            # regex `_(\d+)x(\d+)x(\d+)_r(\d+)` matching. Perminov upstream uses
            # `_m{rank}_` ("m" for the number of multiplications); we rename on
            # download so future re-scans correctly detect this scheme as present.
            out = f"src/main/resources/schemes/section{sec}/perminov-{args.field}_{shape_str}_r{per_field}_{args.field}.json"
            try:
                urllib.request.urlretrieve(url, out)
                print(f"  wrote {out}")
                done += 1
            except Exception as e:
                print(f"  ERROR: {e}", file=sys.stderr)
        else:
            print(f"\n=== {shape_str}: source `{src}` — no handler, skipping ===")
    print(f"\nDone. Imported {done} shapes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
