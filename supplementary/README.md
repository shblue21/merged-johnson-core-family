# Supplementary algebra checks

These scripts accompany the manuscript *Cores of the Merged Johnson Graphs
\(J(2k+1,k)_{\{2,\ldots,k-2\}}\)*. They perform exact symbolic identity
checks and finite-range regression checks.

## Requirements

The scripts were tested with Python 3.9 and SymPy 1.14.0.

    python3 -m venv .venv
    . .venv/bin/activate
    python3 -m pip install -r requirements.txt

## Running the checks

From the paper directory, run:

    make check
    make check-optimized

The optimized-mode target confirms that every check uses explicit failure
handling rather than Python statements that can be disabled with `-O`.

or run the scripts separately:

    python3 supplementary/check_manuscript_identities.py
    python3 supplementary/reconstruct_algebra.py

## Scope

check_manuscript_identities.py checks the following identities used in the
proof symbolically:

- the specialization of the general fibre moment to \(s=2k+3\);
- the even- and odd-module identities for proper divisors;
- the \(D_j=0\) auxiliary identity;
- all four endpoint-interpolation formulas in Appendix A;
- the odd- and even-\(K\) endpoint core identities.

It also performs finite regression checks of the Eberlein endpoint formulas
for \(4\le k\le40\), and derives the parity grids from \(j\) while checking
coefficient nonnegativity and strict \(B_1\)-coefficient positivity for
\(5\le K\le201\).

reconstruct_algebra.py separately solves the defining affine systems and
compares the resulting interpolation coefficients with the factored formulas
in the manuscript. It also performs finite-range checks of:

- Eberlein endpoints and valencies for \(4\le k\le80\);
- candidate fibre sizes for \(4\le k\le500\);
- coefficient signs and strict \(B_1\)-coefficient positivity for
  \(5\le K\le500\);
- the displayed independent-family construction for \(4\le k\le20\);
- endpoint core identities rebuilt from the Eberlein sum for
  \(5\le K\le150\).

The symbolic comparisons are identities in their formal variables. The
finite loops are regression checks only over the ranges stated above. The
scripts do not formalize the extremal upper bound, the maximum-family
classification, the minimum-rank retraction argument, Schur orthogonality,
or the representation-theoretic support reduction.

The code uses X, Y, and Z for the normalized distance-distribution variables
denoted by \(p_1,p_{k-1},p_k\) in the manuscript.

A successful run exits with status zero and ends with:

    All checks in check_manuscript_identities.py completed successfully.
    All reconstruction and finite-range checks completed successfully.
