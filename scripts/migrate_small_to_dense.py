#!/usr/bin/env python3
"""#187: re-serialise SMALL schemes (maxDim < THRESHOLD) from the sparse-map
factor format (u_sparse/v_sparse/w_sparse) back to the dense array format
(u/v/w) for readability. Large schemes (maxDim >= THRESHOLD) keep sparse.

Pure positional densification: dense[k][pos] = sparse_map[k][pos]. The dense
reader (SchemeIO.fromJson) and the sparse reader use identical position indexing
(same m=rank, n=shape, same col-major-W positions), so this round-trips exactly.
All non-factor metadata (fields[], commutative, lineage*, scheduled_additions,
…) is preserved by editing the parsed dict in place. Run with --apply to write;
default is dry-run (verify only)."""
import json, sys, glob, os, collections

THRESHOLD = 16
APPLY = "--apply" in sys.argv

def densify(sparse_map, rank, dim):
    dense = [[0]*dim for _ in range(rank)]
    for kstr, e in sparse_map.items():
        k = int(kstr)
        for pos, coef in zip(e["i"], e["c"]):
            dense[k][pos] = coef
    return dense

def transform(d):
    n = d["n"]; rank = d["m"]
    dimU, dimV, dimW = n[0]*n[1], n[1]*n[2], n[0]*n[2]
    out = collections.OrderedDict()
    for key, val in d.items():
        if key == "u_sparse": out["u"] = densify(val, rank, dimU)
        elif key == "v_sparse": out["v"] = densify(val, rank, dimV)
        elif key == "w_sparse": out["w"] = densify(val, rank, dimW)
        else: out[key] = val
    return out

files = []
for sec in range(2, THRESHOLD):
    files += glob.glob(f"src/main/resources/schemes/section{sec}/*.json")

converted = skipped = 0
for f in sorted(files):
    with open(f) as fh: d = json.load(fh, object_pairs_hook=collections.OrderedDict)
    if "u_sparse" not in d: skipped += 1; continue
    out = transform(d)
    if APPLY:
        with open(f, "w") as fh: json.dump(out, fh, indent=2); fh.write("\n")
    converted += 1

print(f"{'APPLIED' if APPLY else 'DRY-RUN'}: {converted} small files {'converted' if APPLY else 'would convert'} to dense, {skipped} already non-sparse")
