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


def render_win_anatomy():
    """Anatomy of the wins: which derivation mechanism produced each shape where
    our derived rank strictly beats every external holding. Joins the comparison
    JSON with catalog.json's per-scheme lineage_compact."""
    import re
    from collections import Counter

    src = os.path.join(ROOT, 'docs', 'comparison',
                       f'us-vs-fmm-vs-perminov-{COMPARISON_FIELD}.json')
    cat_src = os.path.join(ROOT, 'docs', 'catalog.json')
    comp = json.load(open(src))
    cat = json.load(open(cat_src))['schemes']

    byfmt = {}
    for e in cat:
        byfmt.setdefault(tuple(e['format']), []).append(e)

    def pick(fmt, rank):
        cands = [e for e in byfmt.get(tuple(fmt), [])
                 if e['rank'] == rank and COMPARISON_FIELD in e['fields']
                 and not e.get('commutative')]
        derived = [e for e in cands if e['file'].startswith('derived')]
        return (derived or cands or [None])[0]

    def root_op(lc):
        if lc.startswith('R*['):
            return 'leaf-pair recombination'
        if lc.startswith(('R[', 'Rta[')):
            return 'recombination with allocation'
        depth = 0
        for i, ch in enumerate(lc):
            if ch in '[(':
                depth += 1
            elif ch in '])':
                depth -= 1
            elif depth == 0:
                if lc.startswith(' ⊗ˢ', i):
                    return 'serendipitous product'
                if lc.startswith(' ⊗ ', i):
                    return 'Kronecker product'
                if lc.startswith((' +p ', ' +n '), i):
                    return 'axis concatenation'
                if lc.startswith(' +m ', i):
                    return 'contraction sum'
        if re.search(r'[↓∖]\[', lc):
            return 'downward projection'
        return 'other'

    ALLOC = re.compile(r'^R\*?\[.*?; ([\d,]+) \| ([\d,]+) \| ([\d,]+)[\];]')
    BASE = re.compile(r'^R\*?\[(.+?);')

    roots, bases = Counter(), Counter()
    unbalanced = balanced = 0
    margins = []
    wins = [r for r in comp['rows'] if r.get('_status') == 0]
    for r in wins:
        e = pick(r['format'], r['us_derived'])
        lc = (e or {}).get('lineage_compact')
        if not lc:
            roots['(no lineage)'] += 1
            continue
        roots[root_op(lc)] += 1
        ext = min(x for x in (r.get('fmm'), r.get('perminov'), r.get('us_imported'))
                  if x is not None)
        margins.append(ext - r['us_derived'])
        m = ALLOC.match(lc)
        if m:
            spread = any(max(map(int, g.split(','))) > min(map(int, g.split(',')))
                         for g in m.groups())
            unbalanced += spread
            balanced += not spread
        if lc.startswith(('R[', 'R*[', 'Rta[')):
            b = BASE.match(lc).group(1)
            b = re.sub(r'-[0-9a-f]{7}(\.json)?', '', b)     # strip content hash
            b = re.sub(r'\^ABC->[A-Z]{3}', '', b)           # fold orientation variants
            b = re.sub(r'→⟨[\d,]+⟩', '', b)                 # fold OrientAs wrappers
            bases[b.strip()] += 1

    margins.sort()
    n = len(margins)
    out = [header(["GenerateSourceComparison (Java)", "render_tables.py"], [src, cat_src])]
    out.append("\\begin{table}[h]\n\\centering")
    out.append("\\caption{Anatomy of the " + str(len(wins)) + " strict wins over "
               "$\\mathbb{" + COMPARISON_FIELD + "}$ (Table~\\ref{tab:us-vs-catalogs}): "
               "root operator of the winning scheme's lineage, and the most frequent "
               "combining bases among winning recombinations (orientation variants folded). "
               f"Of the recombination-rooted wins, {unbalanced} use "
               f"\\emph{{unbalanced}} allocations and {balanced} balanced ones. Margins over "
               f"the best external holding: {sum(margins)} ranks in total, median "
               f"{margins[n // 2]}, maximum {margins[-1]}.}}")
    out.append("\\label{tab:win-anatomy}")
    out.append("\\begin{tabular}{lr@{\\qquad}lr}\n\\hline")
    out.append("root operator & wins & top combining base & wins \\\\\n\\hline")
    top_roots = roots.most_common()
    top_bases = bases.most_common(len(top_roots))
    for i in range(max(len(top_roots), len(top_bases))):
        left = f"{top_roots[i][0]} & {top_roots[i][1]}" if i < len(top_roots) else " & "
        if i < len(top_bases):
            bname, bcount = top_bases[i]
            if re.fullmatch(r'[\dx]+', bname):
                bname = "$\\langle " + bname.replace('x', ',') + " \\rangle$"
            right = f"{bname} & {bcount}"
        else:
            right = " & "
        out.append(f"{left} & {right} \\\\")
    out.append("\\hline\n\\end{tabular}\n\\end{table}")
    dst = os.path.join(GEN, 'win-anatomy.tex')
    open(dst, 'w').write("\n".join(out) + "\n")
    print(f"wrote {os.path.relpath(dst, ROOT)} "
          f"({len(wins)} wins; {unbalanced} unbalanced / {balanced} balanced)")


if __name__ == '__main__':
    os.makedirs(GEN, exist_ok=True)
    render_us_vs_catalogs()
    render_win_anatomy()
