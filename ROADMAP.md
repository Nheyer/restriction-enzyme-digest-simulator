# Roadmap

## Continuing an archived project

Upstream (`wl5e/restriction-enzyme-digest-simulator`) was archived on GitHub and
can no longer accept pull requests. This fork continues the project.

Two pull requests were opened against upstream before the archive and cannot be
merged there:

| PR | Work |
|----|------|
| [#2](https://github.com/wl5e/restriction-enzyme-digest-simulator/pull/2) | NEB/REBASE recognition strings, IUPAC ambiguity, non-palindromic and type IIS enzymes |
| [#3](https://github.com/wl5e/restriction-enzyme-digest-simulator/pull/3) | pip-installable, `enzyme-digest` console script |

Both landed here directly instead. The original work is MIT-licensed and the
original author's copyright stands; see [LICENSE](LICENSE).

## Under consideration: a C++ port

Notes toward that decision, so it gets made on evidence rather than instinct.

### What it costs today

Measured on Python 3.12, random sequence, this machine:

| Workload | Time |
|----------|------|
| 6 enzymes over 5 Mbp | 0.56 s (~54 Mbp-enzyme/s) |
| Full 262-enzyme database over 1 Mbp | 5.5 s |
| Full database over a 5 Mbp bacterial genome | ~28 s (extrapolated) |

### Why a straight port probably will not pay

The inner loop is `re.finditer` over a pre-compiled pattern - already C. Porting
the same algorithm to C++ swaps one C regex engine for another, and `std::regex`
is not famously fast. The likely outcome is a rewrite that is no quicker.

The real cost is structural: the tool makes **one full pass per enzyme**, so a
full-database digest scans the sequence 262 times. A single multi-pattern pass -
Aho-Corasick over the expanded sites, or a bit-parallel matcher - collapses that
by roughly the number of enzymes. That is available in Python today, and would
be the honest first move; it also makes any later port a port of a good
algorithm rather than a bad one.

So treat a port as a decision about **deployment** - a single binary, no
interpreter, embedding in a C++ pipeline - and not about speed. If speed is the
goal, fix the pass count first and re-measure.

### If it does go ahead

* `restriction_enzymes.csv` is the data and ports unchanged. Keep it as the
  source of truth; do not re-encode the table as C++ literals, or the two copies
  will drift.
* The cut-coordinate convention is the part most likely to be mis-ported:
  0-based offsets from the start of the site, both strands expressed in the
  top-strand frame, offsets allowed to fall outside the site. A bottom-strand
  hit mirrors as `start + len(site) - p`.
* Port the invariant tests **first**, before any digestion code: fragments
  summing to sequence length, overhang self-complementarity for palindromic
  in-site cutters, and `TestDefiniteMatchingIsPinned`'s 15x15 IUPAC table. They
  are what catch a sign error or an off-by-one in the coordinate frame, and they
  are language-independent.
* `definite` / `possible` ambiguity handling, gap stripping, and origin-spanning
  sites on circular molecules are behaviours the tool promises, not
  implementation details. Port them or document their removal.
* Keep the CLI surface identical so the two implementations can be
  differentially tested against each other on the same FASTA files - the same
  technique used to verify that `--ambiguity definite` preserved the original
  matching rule.
