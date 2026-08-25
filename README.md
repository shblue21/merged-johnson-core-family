# Cores of the Merged Johnson Graphs

This repository contains the manuscript and supplementary algebra checks for

\[
G_k=J(2k+1,k)_{\{2,\ldots,k-2\}},\qquad k\ge4.
\]

The main theorem is

\[
\operatorname{End}(G_k)=\operatorname{Aut}(G_k)\cong S_{2k+1}.
\]

Author: Jihun Kim

Email: jihunkimkw@gmail.com

Repository: https://github.com/shblue21/merged-johnson-core-family

## Contents

- main.tex, references.bib — manuscript source.
- main.pdf — compiled manuscript.
- supplementary/check_manuscript_identities.py — formal-variable identity
  checks and finite regression tests.
- supplementary/reconstruct_algebra.py — separate reconstruction of the
  endpoint-interpolation coefficients and additional finite regression tests.
- supplementary/README.md — exact scope and limitations of the checks.
- requirements.txt — pinned Python dependencies.
- CITATION.cff — citation metadata.
- RELEASE_NOTES.md — notes for the v1.0.0 release candidate.
- LICENSE-CODE — MIT License for software files.
- LICENSE-MANUSCRIPT — CC BY 4.0 for the manuscript and prose documentation.

## Python environment

    python3 -m venv .venv
    . .venv/bin/activate
    python3 -m pip install -r requirements.txt

## Run the algebra checks

    make check
    make check-optimized

The checks complete in a few seconds on a typical laptop. See
supplementary/README.md for the distinction between symbolic identities,
finite-range regression checks, and proof steps not mechanized by the code.

## Build the paper

A TeX installation providing pdflatex and bibtex is required.

    make pdf

To build the PDF and run both normal and optimized-mode checks:

    make verify

For a byte-reproducible release PDF, use a fixed SOURCE_DATE_EPOCH:

    make reproducible SOURCE_DATE_EPOCH=1787616000

## Citation

Citation metadata is provided in CITATION.cff. The Zenodo DOI will be added
after the first archived release.

## License

The Python scripts and related software files are licensed under the MIT
License. The manuscript, bibliography, and prose documentation are licensed
under CC BY 4.0. See LICENSE, LICENSE-CODE, and LICENSE-MANUSCRIPT.
