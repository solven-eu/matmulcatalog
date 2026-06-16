# Paper build + arXiv packaging.
#
#   make pdf      build article.pdf (latexmk: pdflatex + bibtex passes)
#   make arxiv    build the PDF, then bundle an arXiv-ready SOURCE tarball
#                 (article-arxiv.tar.gz). arXiv prefers LaTeX source over a
#                 bare PDF, so this is what you upload.
#   make clean    remove build artifacts and the tarball

PAPER := article
TEX   := $(PAPER).tex
PDF   := $(PAPER).pdf
BBL   := $(PAPER).bbl
ARXIV := $(PAPER)-arxiv.tar.gz

LATEXMK := latexmk -pdf -interaction=nonstopmode -halt-on-error

.PHONY: pdf arxiv arxiv-bundle clean

pdf:
	$(LATEXMK) $(TEX)

# Build the PDF (which generates the .bbl), then bundle the source.
arxiv: pdf arxiv-bundle

# Bundle every LaTeX source under paper/ (recursively: macros, sections/,
# tables/, generated/, theory/, …) plus refs.bib and the generated .bbl, so a
# fresh checkout on arXiv resolves all \input{} and \cite{} without needing the
# local working tree. Build artifacts (.aux/.log/.pdf/.fls/…) are excluded.
# Standalone (no PDF prerequisite) so CI can reuse the .bbl a prior compile
# step already produced instead of running latexmk twice.
arxiv-bundle:
	@echo "Bundling arXiv source tarball from $(PAPER).tex …"
	@tar -czf $(ARXIV) \
	  $$(find . -type f \( -name '*.tex' -o -name '*.bib' -o -name '*.bbl' -o -name '*.cls' -o -name '*.sty' \) \
	     | sed 's|^\./||' | sort)
	@echo "Wrote $(ARXIV) ($$(du -h $(ARXIV) | cut -f1)); contents:"
	@tar -tzf $(ARXIV) | sed 's/^/    /'

clean:
	latexmk -C $(TEX)
	rm -f $(ARXIV)
