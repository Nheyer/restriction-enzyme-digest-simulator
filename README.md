# Restriction Enzyme Digest Simulator

A CLI tool to simulate restriction enzyme digests on DNA sequences, compute fragment sizes, and generate an ASCII gel electrophoresis visualization.

## Features
- In silico digestion of linear or circular DNA from FASTA files.
- Supports multiple enzymes simultaneously (combining their cuts).
- Built-in library of 88 restriction enzymes, with specificities taken from [REBASE](http://rebase.neb.com).
- Enzymes are defined in standard NEB/REBASE notation: `G^AATTC`, `GGTCTC(1/5)`.
- IUPAC ambiguity codes in recognition sites (`GTYRAC`, `CCANNNNNNTGG`).
- Non-palindromic enzymes: both strands are searched, and bottom-strand sites cut where the enzyme actually cuts them.
- Type IIS enzymes that cut outside their recognition site (BsaI, BbsI, FokI, ...).
- Filter fragments below a minimum size.
- Output as a table of fragment lengths or ASCII gel with a 100 bp DNA ladder.
- Handles overlapping recognition sites and no-cut scenarios.
- Works on multi-sequence FASTA files, processing each entry separately.

## Installation

Requires Python 3.8 or later. No runtime dependencies.

```bash
pip install git+https://github.com/gorgerat/restriction-enzyme-digest-simulator.git
```

That puts an `enzyme-digest` command on your PATH:

```bash
enzyme-digest --fasta plasmid.fasta --enzymes EcoRI,BamHI
```

From a clone, for development:

```bash
git clone https://github.com/gorgerat/restriction-enzyme-digest-simulator.git
cd restriction-enzyme-digest-simulator
pip install -e ".[test]"     # editable, with pytest
```

Installing is optional - the script runs straight from a checkout with no
install step and no dependencies:

```bash
python enzyme_digest.py --fasta plasmid.fasta --enzymes EcoRI,BamHI
```

The two invocations are equivalent; the examples below use `python enzyme_digest.py`.

### As a library

```python
from enzyme_digest import ENZYME_DB, find_sites, digest_linear, normalise_cuts

enzyme = ENZYME_DB["BsaI"]
hits = find_sites(sequence, enzyme)
cuts, uncut = normalise_cuts([h.cut for h in hits], len(sequence), circular=False)
fragments = digest_linear(len(sequence), cuts, min_fragment=1)
```

## Project layout

| File | Contents |
|------|----------|
| `enzyme_digest.py` | The simulator: parsing, site finding, digestion, output, CLI |
| `pyproject.toml` | Packaging metadata and the `enzyme-digest` entry point |
| `enzyme_data.py` | Static reference data only - IUPAC codes, the complement table, the enzyme table, the ladder |
| `tests/` | pytest suite |

`enzyme_data.py` imports nothing from the simulator, so the enzyme table can be regenerated or diffed on its own.

## Usage

```bash
python enzyme_digest.py --fasta <file.fasta> --enzymes <enzyme1,enzyme2,...> [options]
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-f`, `--fasta` | Path to input FASTA file | (required) |
| `-e`, `--enzymes` | Comma-separated enzyme names or custom specs | (required) |
| `-c`, `--circular` | Treat DNA as circular | `False` (linear) |
| `-m`, `--min-fragment` | Minimum fragment length to report (bp) | 1 |
| `-o`, `--output` | Output mode: `table`, `gel`, or `both` | `table` |
| `--gel-height` | Height of the ASCII gel | 30 |
| `--list-enzymes` | List the built-in enzymes and exit | |
| `--convert SPEC` | Show an enzyme in both notations and exit | |
| `--version` | Print the version and exit | |

### Enzyme specifications

Built-in enzymes are given by name (case-insensitive): `EcoRI`, `BamHI`, `BsaI`, `HincII`, ... Run `--list-enzymes` for the full list with recognition sites and end types.

Custom enzymes can be written in any of these forms:

| Form | Example | Meaning |
|------|---------|---------|
| `NAME:RECOGNITION` | `PmeI:GTTT^AAAC` | NEB caret notation - the cut is marked inside the site |
| `NAME:RECOGNITION` | `MyI:GGTCTC(1/5)` | NEB offset notation - cuts 1 nt (top) and 5 nt (bottom) past the site |
| `NAME:SEQ:OFFSET` | `PmeI:GTTTAAAC:4` | Top-strand cut offset; the bottom strand is assumed to mirror it |
| `NAME:SEQ:TOP:BOTTOM` | `MyI:GGTCTC:7:11` | Both cut offsets, measured from the start of the site |

A bare recognition string with no name also works (`-e "G^AATTC"`).

Recognition sites may use the IUPAC ambiguity codes `RYSWKMBDHVN`. A code in the site matches a base in the target only when that base is unambiguously one the site accepts, so an `N` in the enzyme matches an `N` in the sequence, but an `A` does not.

### Converting between notations

```bash
$ python enzyme_digest.py --convert BsaI
Enzyme:            BsaI
Recognition site:  GGTCTC(1/5)
NAME:SEQ:OFFSET:   not representable - BsaI cuts outside its recognition site ...
Cut (top/bottom):  7/11 from the start of the site
Ends:              5' overhang, 4 nt
Palindromic:       no

$ python enzyme_digest.py --convert HincII:GTYRAC:3
Enzyme:            HincII
Recognition site:  GTY^RAC
NAME:SEQ:OFFSET:   HincII:GTYRAC:3
Cut (top/bottom):  3/3 from the start of the site
Ends:              blunt
Palindromic:       yes
Ambiguity codes:   4 concrete sites (GTCAAC, GTCGAC, GTTAAC, GTTGAC)
```

`NAME:SEQ:OFFSET` records only the top-strand cut and assumes the bottom-strand cut mirrors it, so it cannot describe an enzyme that cuts asymmetrically or outside its site. The converter says so instead of quietly dropping the bottom-strand cut.

### Examples

1. **Linear digest with two enzymes, table output**
   ```bash
   python enzyme_digest.py -f plasmid.fasta -e EcoRI,BamHI
   ```

2. **Circular digest, gel visualization**
   ```bash
   python enzyme_digest.py -f circular.fasta -e EcoRI,BamHI --circular --output gel
   ```

3. **Single custom enzyme and filtering small fragments**
   ```bash
   python enzyme_digest.py -f dna.fasta -e PmeI:GTTTAAAC:4 --min-fragment 50
   ```

4. **Run with multiple sequences in one FASTA**
   ```bash
   python enzyme_digest.py -f multi.fasta -e HindIII,NotI --output both
   ```

5. **Golden Gate enzyme on both strands**
   ```bash
   python enzyme_digest.py -f construct.fasta -e BsaI
   ```
   Sites are reported per strand, since BsaI is non-palindromic and cuts 1/5 nt downstream of its site.

## Cut position conventions

Cut positions are 0-based offsets from the start of the recognition site, measured along the top strand: the cut falls immediately before the base at that offset.

- `cut_top` is where the top strand is cut, `cut_bottom` where the bottom strand is cut, both in the same top-strand frame.
- `cut_bottom - cut_top` is the overhang: positive for a 5' overhang, 0 for a blunt cut, negative for a 3' overhang.
- Both offsets can lie outside the site. Type IIS enzymes cut downstream of it, and `^GATC` (MboI) cuts immediately before it.
- A site found on the bottom strand is mirrored: an offset `p` in the enzyme's own frame lands at `start + len(site) - p` in the sequence.
- On linear DNA a cut that falls beyond either end is reported as a recognised-but-uncut site. On circular DNA it wraps around.

## Known limitations

- Enzymes that cut on both sides of their recognition site (BaeI, BcgI, CspCI, `(10/15)ACNNNNGTAYC(12/7)`) are rejected rather than partly parsed: two cuts per site do not fit the single-cut fragment model.
- Recognition sites that span the origin of a circular sequence are not found; only cut positions wrap.
- Fragment lengths are top-strand lengths, so overhangs are not reflected in the reported sizes.
- Methylation sensitivity, star activity, and enzymes requiring two sites are not modelled.

## Testing

```bash
pip install -e ".[test]"
pytest tests/ -v
```

Or without installing anything:

```bash
pip install pytest
PYTHONPATH=. pytest tests/ -v
```

## License

MIT License – see [LICENSE](LICENSE).

---
*Author: Collins Amatu Gorgerat, 2026*
