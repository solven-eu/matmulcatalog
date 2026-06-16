# GitHub Pages source

This directory is the publish root for the
[strassen catalog browser](https://solven-eu.github.io/matmulcatalog/) — a
static HTML page that lets you filter the ~2,300-scheme catalog by
field, max-dimension, source, and shape.

## Files

| file | role |
|---|---|
| `index.html` | UI: filter controls + results table |
| `style.css` | basic styling, single column on mobile |
| `catalog.js` | vanilla-JS loader + filter logic |
| `catalog.json` | auto-generated manifest, **regenerate after scheme changes** |

## Regenerating `catalog.json`

```bash
mvn -q test-compile
CP=$(mvn dependency:build-classpath -q -Dmdep.outputFile=/dev/stdout 2>/dev/null)
java -cp target/classes:target/test-classes:$CP \
     io.cormoran.strassen.v3.catalog.GenerateCatalogManifest
```

That walks every `*.json` under `src/main/resources/schemes/section{N}/`,
reads each via `SchemeIO`, computes the addition count, and writes a
flat JSON list of `{format, max_dim, field, rank, additions, source,
verified, complex, file}`.

## Enabling on GitHub Pages

In the repo settings → Pages → set source to *Deploy from a branch* with
**branch=master** and **folder=/docs**. The site becomes available at
`https://solven-eu.github.io/matmulcatalog/`.

## Local preview

Just `python3 -m http.server -d docs 8000` and visit
`http://localhost:8000/`.
