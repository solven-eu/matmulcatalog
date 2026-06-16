# Paper draft

LaTeX source for the paper version of the catalog. Roadmap entry:
`ROADMAP.md` → "Publishable article in LaTeX".

## Files

```
paper/
├── article.tex                # main document, \input's sections
├── macros.tex                 # shape-notation macros (\nmpfield, etc.)
├── refs.bib                   # bibliography
├── sections/
│   ├── abstract.tex
│   ├── intro.tex
│   ├── notation.tex           # field discipline & shape conventions
│   ├── architecture.tex       # catalog architecture (explicit / cited / derived)
│   ├── strategies.tex         # composition strategies — the meaty section
│   ├── nonoverlap.tex         # non-overlap property (Theorem)
│   ├── search.tex             # frontier-closure search algorithm
│   ├── tables.tex             # comparison tables scaffold
│   └── openquestions.tex      # what we can't answer yet
└── tables/
    ├── nc-Q.tex   nc-R.tex   nc-C.tex   nc-F2.tex
    └── c-Q.tex    c-R.tex
```

## Build

```
cd paper
latexmk -pdf article.tex          # convenient if you have latexmk
# or:
pdflatex article && bibtex article && pdflatex article && pdflatex article
```

Tested against TeX Live 2024. No exotic packages -- everything from
the standard texlive-latex-extra set.

## SOTA conventions

Several edge cases require deterministic decisions before the
comparison tables auto-regenerate. Currently captured in
`paper/sota-conventions.md` (TODO — write this) — at minimum:

- mixed-field (Q-coefficient scheme tagged for R)
- border-rank vs exact rank (border has its own appendix)
- partially verified cited bounds (dagger)
- F2 promotion from Z (not yet automated)

## Auto-regen tooling

Tracked as a separate roadmap item. The intent is a small Java
driver that reads `docs/catalog.json`, `docs/cited-bounds.json`,
`docs/derived-from-cited-bounds.json` and emits `paper/tables/*.tex` files
ready to be `\input`'d.

Today's `paper/tables/*.tex` files are populated by hand and will
rot on every catalog change; treat the values as illustrative until
the regen pipeline lands.

## Editing conventions

- Always use the shape-notation macros from `macros.tex`. Never
  write `<n,m,p>` or `(n,m,p)` directly — pick `\nmpfield{R}{n}{m}{p}`
  or its commutative / multiplications-count / addition-count
  variants.
- One paragraph = one idea. Section files should remain small enough
  that the section table-of-contents reads as the section's outline.
- Don't auto-commit when a value in a table is uncertain; flag with
  `\dag` and a footnote.
- New citations go into `refs.bib` with a key matching the
  `author{year}` pattern (`strassen1969`, `alphatensor2022`, etc.).
