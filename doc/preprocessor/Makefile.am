if HAVE_PDFLATEX
if HAVE_BEAMER
pdf-local: preprocessor.pdf
endif
endif

SRC = preprocessor.tex expr.png expr-sharing.png matrices.png overview.png json-preprocessor.png

EXTRA_DIST = $(SRC)

preprocessor.pdf: $(SRC)
	$(PDFLATEX) preprocessor
	$(PDFLATEX) preprocessor

clean-local:
	rm -f *.pdf *.toc *.aux *.log *.nav *.snm *.vrb *.out *~
