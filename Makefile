PYTHON ?= python3
PDFLATEX ?= pdflatex
BIBTEX ?= bibtex
SOURCE_DATE_EPOCH ?= 1787616000

.PHONY: all pdf check check-optimized verify reproducible clean

all: main.pdf

pdf: main.pdf

main.pdf: main.tex references.bib
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error main.tex
	$(BIBTEX) main
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error main.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error main.tex

check:
	$(PYTHON) supplementary/check_manuscript_identities.py
	$(PYTHON) supplementary/reconstruct_algebra.py

check-optimized:
	$(PYTHON) -O supplementary/check_manuscript_identities.py
	$(PYTHON) -O supplementary/reconstruct_algebra.py

verify:
	$(MAKE) clean
	$(MAKE) all
	$(MAKE) check
	$(MAKE) check-optimized

reproducible:
	$(MAKE) clean
	SOURCE_DATE_EPOCH=$(SOURCE_DATE_EPOCH) $(MAKE) all

clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.pdf
