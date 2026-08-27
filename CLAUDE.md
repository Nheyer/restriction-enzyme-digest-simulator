# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

A CLI tool for in silico restriction digests. This repo is a **fork continuing an
archived project** — upstream (`wl5e/restriction-enzyme-digest-simulator`) is
archived and fully read-only, so its PRs cannot be commented on, closed, or
merged. All work lands here on `main`.

A **C++ port is under consideration**. [ROADMAP.md](ROADMAP.md) has the
benchmarks behind that decision and what a port must preserve — read it before
starting one. Short version: the inner loop is already C, and the real cost is
one scan pass per enzyme, so the win is algorithmic rather than linguistic.

## Commands

```bash
pip install -e ".[test]"          # editable install, with pytest
pytest tests/ -v                  # whole suite (134 tests)
pytest tests/ -k TestFindSites    # one class
pytest tests/test_enzyme_digest.py::TestFindSites::test_two_sites   # one test

python tools/build_enzyme_table.py            # regenerate tables from the CSV
python tools/build_enzyme_table.py --check    # CI gate: fails if out of sync

python enzyme_digest.py -f seq.fasta -e EcoRI,BsaI   # from a checkout
enzyme-digest -f seq.fasta -e EcoRI,BsaI             # installed console script
```

CI runs on 3.8 / 3.9 / 3.12, and smoke-tests the console script from outside the
source tree — running it from the repo root would exercise the checkout rather
than the installed copy.

## Architecture

**Two flat modules, deliberately not a package.** `enzyme_digest.py` is all the
logic and the CLI; `enzyme_data.py` is lookup tables only. Both are declared in
`pyproject.toml` under `py-modules`, so `python enzyme_digest.py` keeps working
from a checkout while `pip install` still ships both. A test asserts
`enzyme_data` defines no functions or classes and does not import
`enzyme_digest`; keep that direction of dependency.

**`restriction_enzymes.csv` is the source of truth for enzyme data.** Edit the
CSV, then regenerate — never hand-edit the tables in `enzyme_data.py`, which sit
between `BEGIN/END GENERATED TABLES` markers and are overwritten. The CSV uses
NEB bench notation, `^` for the top-strand cut and `_` for the bottom:

```
EcoRI,G^AATTC        AatII,G_ACGT^C        BsaI,GGTCTCN^NNNN_
```

Generation (rather than reading the CSV at import) exists for two reasons:
`enzyme_data` must stay pure data, and setuptools ships `py-modules` but not
data files beside them, so a runtime CSV read would not survive `pip install`.
`--enzyme-db` loads a different CSV at runtime.

**Cut coordinates are the thing to get right.** Both offsets are 0-based from
the first base of the recognition site, both measured along the *top* strand,
with the cut falling immediately before that offset. `cut_bottom - cut_top` is
the overhang: positive 5', zero blunt, negative 3'. Offsets may fall outside the
site — type IIS cuts downstream, `^GATC` (MboI) cuts before. A site found on the
bottom strand mirrors: enzyme-frame offset `p` lands at `start + len(site) - p`,
so a reverse-orientation hit cuts *upstream* of where it was found.

Cuts landing outside the sequence are dropped on linear DNA (reported as
recognised-but-uncut) and wrapped modulo length on circular. Under `--circular`
the scan window is extended by `len(site) - 1` so origin-spanning sites are
found; no site can be double-counted because the appended tail is shorter than
the site.

**Ambiguity has two modes and an input-alphabet policy**, all in
`site_to_regex`. `definite` matches only when the target base set is a subset of
what the site accepts ("where will this definitely cut?"); `possible` needs only
an intersection ("could this cut at all?"). Matching uses a regex built from the
site — not string expansion, which explodes on sites like XcmI's nine Ns.
`expand_ambiguity()` exists for callers that want the enumeration. Without
`--ambiguity`, only ACGT is accepted and anything else is a hard error, because
a silently unmatched ambiguous base is how a sequence gets wrongly cleared as
uncut.

**Enzyme classes the model cannot represent are quarantined, not dropped.**
Nicking enzymes (`Nb.*`/`Nt.*`) and type IIB dual-cut enzymes go into separate
tables so they produce a named error rather than "unknown enzyme". A single cut
coordinate per site is baked into the fragment model.

## Conventions that matter here

**Enzyme specificities come from REBASE, never from recall.** Wrong cut data is
silent and biological — the worst failure mode in this repo. NEB's web list
returns 403 to fetches; use `http://rebase.neb.com/rebase/link_emboss_e`
(`name pattern len ncuts blunt c1 c2 c3 c4`, where a 0-based offset is
`c if c > 0 else c + 1`).

**Do not derive test expectations from the implementation.** Several tests exist
specifically to be independent of the code they check, and re-deriving them with
the same comprehension would make them test themselves:

- `TestDefiniteMatchingIsPinned` — a hand-written 15-code IUPAC table.
- overhang self-complementarity for palindromic in-site cutters — forced by
  symmetry, so it is independent of REBASE too.
- fragments summing to sequence length across every enzyme, linear and circular
  — this is what catches sign errors and off-by-ones in bulk.

When changing matching or coordinate code, a differential test against the
previous revision (same corpus, both versions, compare hits) is the strongest
check available; include `B/D/H/V` in any ambiguity corpus, since those are
where subset and intersection diverge most.

## Repo facts

`origin` must stay on **HTTPS** — this machine has no SSH key, and `gh auth
setup-git` installs a credential helper for `https://github.com` only, which an
SSH remote never consults.
