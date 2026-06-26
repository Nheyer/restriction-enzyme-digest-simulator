# Restriction Enzyme Digest Simulator

A CLI tool to simulate restriction enzyme digests on DNA sequences, compute fragment sizes, and generate an ASCII gel electrophoresis visualization.

## Features
- In silico digestion of linear or circular DNA from FASTA files.
- Supports multiple enzymes simultaneously (combining their cuts).
- Built-in library of common restriction enzymes.
- Custom enzyme definition: `NAME:RECOGNITION_SEQUENCE:CUT_OFFSET`.
- Filter fragments below a minimum size.
- Output as a table of fragment lengths or ASCII gel with a 100 bp DNA ladder.
- Handles overlapping recognition sites and no-cut scenarios.
- Works on multi-sequence FASTA files, processing each entry separately.

## Installation

Requires Python 3.8 or later.

```bash
git clone https://github.com/gorgerat/restriction-enzyme-digest-simulator.git
cd restriction-enzyme-digest-simulator
```

No external runtime dependencies. Install `pytest` for running tests.

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

### Enzyme specifications

- **Built-in**: `EcoRI`, `BamHI`, `HindIII`, `EcoRV`, `PstI`, `SmaI`, `XbaI`, `NotI`, `SacI`, `KpnI`, `SpeI`, `BglII`, `NcoI`, `NdeI`, `XhoI`, `SalI`, `ClaI`, `HpaI`, `MluI`, `NheI`
- **Custom**: `NAME:SEQ:OFFSET` – e.g., `PmeI:GTTTAAAC:4` (cuts between T and A, 4 bases into the site)

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

## Testing

```bash
pip install pytest
PYTHONPATH=. pytest tests/ -v
```

## License

MIT License – see [LICENSE](LICENSE).

---
*Author: Collins Amatu Gorgerat, 2026*
