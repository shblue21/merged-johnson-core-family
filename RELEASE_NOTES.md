# Release notes

## v1.0.0 — 25 August 2026

Initial public release candidate for the manuscript *Cores of the Merged
Johnson Graphs \(J(2k+1,k)_{\{2,\ldots,k-2\}}\)* and its supplementary
algebra checks.

### Mathematical result

For every integer \(k\ge4\), the manuscript proves

\[
\operatorname{End}(G_k)=\operatorname{Aut}(G_k)\cong S_{2k+1}.
\]

It also determines \(\alpha(G_k)=2k+3\), classifies the maximum independent
sets, and derives the corresponding odd-\(n\) qualitative-independence
hypergraph corollary.

### Supplementary checks

- formal-variable checks of the load-bearing rational identities;
- separate reconstruction of all endpoint-interpolation coefficients;
- explicit normal- and optimized-mode failure handling;
- finite-range Eberlein, divisor, coefficient-sign, construction, and
  endpoint-identity regression tests;
- pinned Python dependencies and a Python 3.9/3.12 CI matrix;
- a deterministic PDF build target based on SOURCE_DATE_EPOCH.

The exact scope and limitations of the programs are documented in
supplementary/README.md.

The deterministic release-candidate PDF has SHA-256
`c00f1a7e2dd16605ef6d52df00583cb89423e96f1ff6195fe8ef090b0a24ad42`.

### Archival identifiers

- Version DOI: https://doi.org/10.5281/zenodo.22096460
- All-versions DOI: https://doi.org/10.5281/zenodo.22096459

### Reproduction

    python3 -m venv .venv
    . .venv/bin/activate
    python3 -m pip install -r requirements.txt
    make verify
    make reproducible

### Licenses

- Software: MIT License.
- Manuscript and prose documentation: CC BY 4.0.
