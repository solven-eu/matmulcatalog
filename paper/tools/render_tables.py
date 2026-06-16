#!/usr/bin/env python3
"""RENDER layer of the paper-table pipeline: docs/comparison/*.json → paper/generated/*.tex.

Pipeline contract (three decoupled layers — see GenerateSourceComparison):
  [1] Java generators           → docs/comparison/*.json   (data; LaTeX-unaware)
  [2] THIS script               → paper/generated/*.tex     (fragments; machine-owned)
  [3] paper/sections/*.tex      → \\input{generated/...}    (hand-written prose)

paper/generated/ is MACHINE TERRITORY: every emitted file carries a DO-NOT-EDIT
header naming the generator chain, the input file and its content hash — so a
reader of the .tex can trace any number back to the code that produced it, and
staleness is detectable by re-hashing.

Run:  python3 paper/tools/render_tables.py     (from the repo root)
"""
import hashlib
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GEN = os.path.join(ROOT, 'paper', 'generated')


def sha7(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()[:7]


def header(generator_chain, inputs):
    lines = ["% AUTO-GENERATED — DO NOT EDIT (machine-owned: paper/generated/)."]
    lines.append("% Pipeline: " + " -> ".join(generator_chain))
    for p in inputs:
        lines.append(f"% Input: {os.path.relpath(p, ROOT)}  (sha256:{sha7(p)})")
    lines.append("% Regenerate: mvn exec:java -Dexec.mainClass=eu.solven.matmul.docs.generate.GenerateSourceComparison"
                 " && python3 paper/tools/render_tables.py")
    return "\n".join(lines) + "\n"


# The comparison is FIELD-SCOPED (one JSON per field — see GenerateSourceComparison).
# The paper renders the Q table: \mathbb{Q} is the field where both external
# catalogs (FMM-Lille, Perminov) carry the most data, so the head-to-head is
# meaningful. 'this work' is our OWN derived best (us_derived); us_imported shows
# whether we also carry the known scheme.
COMPARISON_FIELD = 'Q'


def render_us_vs_catalogs():
    src = os.path.join(ROOT, 'docs', 'comparison',
                       f'us-vs-fmm-vs-perminov-{COMPARISON_FIELD}.json')
    data = json.load(open(src))
    # _status: 0=win (our DERIVED work strictly best over the field), 1=loss,
    # 2=tie, 3=carried. The paper table shows our wins.
    rows = [r for r in data['rows'] if r.get('_status') == 0]

    # Presentation choice (renderer-owned): order by margin over the best
    # external bound, cap the table, keep it readable.
    def ext_best(r):
        vals = [x for x in (r.get('fmm'), r.get('perminov')) if x is not None]
        return min(vals) if vals else None

    def margin(r):
        ext = ext_best(r)
        return (ext - r['us_derived']) if ext is not None else 10**9
    rows.sort(key=margin, reverse=True)
    cap = 24
    shown = rows[:cap]

    out = [header(["GenerateSourceComparison (Java)", "render_tables.py"], [src])]
    out.append("\\begin{table}[h]\n\\centering")
    out.append("\\caption{Shapes where this work's own \\emph{derived} non-commutative rank "
               "over $\\mathbb{" + COMPARISON_FIELD + "}$ strictly improves on every external "
               "holding (FMM-Lille, Perminov). The \\emph{carry} column is the best rank we hold "
               "as an import over the same field. "
               f"Top {len(shown)} of {len(rows)} such shapes by margin; full data in "
               "\\code{docs/comparison/us-vs-fmm-vs-perminov-" + COMPARISON_FIELD + ".json}.}")
    out.append("\\label{tab:us-vs-catalogs}")
    out.append("\\begin{tabular}{lrrrr}\n\\hline")
    out.append("shape & this work & carry & FMM-Lille & Perminov \\\\\n\\hline")
    for r in shown:
        f = r['format']
        def cell(v):
            return str(v) if v is not None else '---'
        out.append(f"$\\langle {f[0]},{f[1]},{f[2]} \\rangle$ & "
                   f"\\textbf{{{r['us_derived']}}} & {cell(r.get('us_imported'))} & "
                   f"{cell(r.get('fmm'))} & {cell(r.get('perminov'))} \\\\")
    out.append("\\hline\n\\end{tabular}\n\\end{table}")
    dst = os.path.join(GEN, 'us-vs-catalogs.tex')
    open(dst, 'w').write("\n".join(out) + "\n")
    print(f"wrote {os.path.relpath(dst, ROOT)} ({len(shown)} rows shown / {len(rows)} wins over {COMPARISON_FIELD})")


if __name__ == '__main__':
    os.makedirs(GEN, exist_ok=True)
    render_us_vs_catalogs()
