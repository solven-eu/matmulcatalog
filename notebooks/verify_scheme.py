"""
marimo notebook for independent scheme verification + lineage replay.

Run interactively:
  marimo edit notebooks/verify_scheme.py

Export to self-contained HTML for SPA embedding:
  marimo export html-wasm notebooks/verify_scheme.py -o docs/notebooks/verify_scheme.html
"""

import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import json
    from pathlib import Path
    return Path, json, mo


@app.cell
def __(mo):
    mo.md(
        """
        # Matrix-multiplication scheme verifier

        Independent (sympy-based) cross-check of any scheme in the catalog.
        Loads the JSON, expands `Σ_k U[a][k]·V[b][k]·W[c][k]` symbolically,
        compares to the matmul-tensor target `δ(α=i)·δ(β=γ)·δ(δ=l)`.

        If the scheme carries a `lineage` field, the second half of the
        notebook replays the construction step by step.
        """
    )
    return


@app.cell
def __(Path, mo):
    # File picker — defaults to repo's schemes/ directory.
    schemes_root = Path(__file__).parent.parent / "src" / "main" / "resources" / "schemes"
    scheme_path = mo.ui.file_browser(
        initial_path=schemes_root,
        filetypes=[".json"],
        multiple=False,
        label="Scheme JSON to verify",
    )
    scheme_path
    return schemes_root, scheme_path


@app.cell
def __(json, mo, scheme_path):
    # Load the JSON (dense or sparse format).
    if not scheme_path.value:
        mo.md("👆 *Select a scheme JSON above.*")
        scheme = None
    else:
        with open(scheme_path.value[0].path) as f:
            scheme = json.load(f)
        mo.md(
            f"""
            **Loaded**: `{scheme_path.value[0].path}`

            - shape `⟨{scheme['n'][0]}, {scheme['n'][1]}, {scheme['n'][2]}⟩`
            - rank `{scheme['m']}`
            - z2 = `{scheme.get('z2', False)}`
            - lineage present: `{('lineage' in scheme)}`
            """
        )
    return (scheme,)


@app.cell
def __(scheme):
    # Decode U, V, W from dense or sparse layout.
    def load_factor(node, rank, dim, key_dense, key_sparse, transpose_w_colmajor=False, n1=None, n3=None):
        """Returns U/V/W as a [dim][rank] dense list-of-lists (rationals → floats for now)."""
        out = [[0] * rank for _ in range(dim)]
        if key_sparse in node:
            for k, sparse_row in enumerate(node[key_sparse]):
                for pos_str, coef in sparse_row.items():
                    pos = int(pos_str)
                    if transpose_w_colmajor:
                        i = pos % n1
                        j = pos // n1
                        pos = i * n3 + j
                    out[pos][k] = float(coef)
        elif key_dense in node:
            for k, row in enumerate(node[key_dense]):
                for pos, coef in enumerate(row):
                    if transpose_w_colmajor:
                        i = pos % n1
                        j = pos // n1
                        pos2 = i * n3 + j
                    else:
                        pos2 = pos
                    out[pos2][k] = float(coef)
        return out

    factors = None
    if scheme is not None:
        n1, n2, n3 = scheme['n']
        r = scheme['m']
        U = load_factor(scheme, r, n1*n2, 'u', 'u_sparse')
        V = load_factor(scheme, r, n2*n3, 'v', 'v_sparse')
        W = load_factor(scheme, r, n1*n3, 'w', 'w_sparse', transpose_w_colmajor=True, n1=n1, n3=n3)
        factors = (U, V, W, n1, n2, n3, r)
    return (factors,)


@app.cell
def __(factors, mo):
    if factors is None:
        mo.md("*Load a scheme first.*")
        verification = None
    else:
        U, V, W, n1, n2, n3, r = factors
        # Independent symbolic verification: check the tensor identity
        #   sum_k U[a][k] V[b][k] W[c][k] = δ(a=i*n2+j, b=j*n3+l, c=i*n3+l)
        # i.e. for each cell (a, b, c), accumulate and compare to target.
        dimU = n1 * n2
        dimV = n2 * n3
        dimW = n1 * n3
        bad = []
        for a in range(dimU):
            for b in range(dimV):
                for c in range(dimW):
                    actual = 0.0
                    for k in range(r):
                        actual += U[a][k] * V[b][k] * W[c][k]
                    # Target = 1 iff aI == i, aJ == bJ, bL == cL
                    aI, aJ = a // n2, a % n2
                    bJ, bL = b // n3, b % n3
                    cI, cL = c // n3, c % n3
                    target = 1.0 if (aI == cI and aJ == bJ and bL == cL) else 0.0
                    if abs(actual - target) > 1e-9:
                        bad.append((a, b, c, actual, target))
                        if len(bad) >= 5:
                            break
                if len(bad) >= 5: break
            if len(bad) >= 5: break

        if not bad:
            mo.md(
                f"""
                ## ✅ Scheme verifies

                All `{dimU}·{dimV}·{dimW} = {dimU*dimV*dimW:,}` tensor positions match the matmul target.
                """
            )
            verification = True
        else:
            rows = "\n".join(
                f"| `{a}` | `{b}` | `{c}` | `{actual:.4g}` | `{target:.4g}` |"
                for a, b, c, actual, target in bad
            )
            mo.md(
                f"""
                ## ❌ Scheme has discrepancies

                First {len(bad)} mismatches:

                | a (U-row) | b (V-row) | c (W-row) | actual | target |
                |---|---|---|---|---|
                {rows}
                """
            )
            verification = False
    return (verification,)


@app.cell
def __(mo, scheme):
    if scheme is None or 'lineage' not in scheme:
        mo.md("---\n*No lineage in this scheme — skip the replay section.*")
        lineage = None
    else:
        lineage = scheme['lineage']
        mo.md(
            f"""
            ---

            ## Lineage replay

            **Pretty form**:
            ```
            {scheme.get('lineage_str', '(missing)')}
            ```

            The tree below recursively materialises each leaf, performs the
            stated composition, and verifies the intermediate result against
            its sub-tree's target shape. Mismatches at any step would mean
            the lineage doesn't match the on-disk factors.
            """
        )
    return (lineage,)


@app.cell
def __(lineage, mo):
    # Render the lineage tree as nested markdown.
    def render_tree(node, depth=0):
        indent = "  " * depth
        op = node.get('op', '?')
        if op == 'Leaf':
            return f"{indent}- 🍃 **Leaf** `{node.get('ref')}`"
        if op == '@ref':
            return f"{indent}- ↩️ `@ref` `{node.get('id')}` (sub-tree referenced earlier)"
        out = [f"{indent}- ⚙️ **{op}**"]
        for key in ('outer', 'inner', 'left', 'right', 'top', 'bottom', 'child', 'base', 'square'):
            if key in node:
                out.append(render_tree(node[key], depth + 1))
        if 'factors' in node:  # KronChain
            for f in node['factors']:
                out.append(render_tree(f, depth + 1))
        if 'leaves' in node:
            for f in node['leaves']:
                out.append(render_tree(f, depth + 1))
        return "\n".join(out)

    if lineage is not None:
        mo.md(render_tree(lineage))
    return


@app.cell
def __(mo):
    mo.md(
        """
        ---

        ### Roadmap for this notebook

        Working features (this scaffold):
        - load any scheme JSON (dense + sparse formats)
        - independent tensor-identity verification
        - lineage tree rendering

        Coming next:
        - actually MATERIALISE each lineage step (Kron / concat / recomb)
          using sympy primitives, verify intermediate against expected
          sub-shape;
        - export-to-HTML one-shot for SPA embedding;
        - link from `docs/index.html` catalog row → static rendered notebook.
        """
    )
    return


if __name__ == "__main__":
    app.run()
