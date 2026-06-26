#!/usr/bin/env python3
"""Restriction Enzyme Digest Simulator."""

import argparse
import sys
import math
import re
from typing import List, Tuple, Dict, Optional

# Common restriction enzymes (name, recognition_seq, cut_offset_0based)
# Cut offset indicates the last nucleotide before the break on the top strand.
ENZYME_DB = {
    "EcoRI": ("GAATTC", 1),   # G^AATTC
    "BamHI": ("GGATCC", 1),   # G^GATCC
    "HindIII": ("AAGCTT", 1), # A^AGCTT
    "EcoRV": ("GATATC", 3),   # GAT^ATC
    "PstI": ("CTGCAG", 5),    # CTGCA^G
    "SmaI": ("CCCGGG", 3),    # CCC^GGG
    "XbaI": ("TCTAGA", 1),    # T^CTAGA
    "NotI": ("GCGGCCGC", 2),  # GC^GGCCGC
    "SacI": ("GAGCTC", 5),    # GAGCT^C
    "KpnI": ("GGTACC", 5),    # GGTAC^C
    "SpeI": ("ACTAGT", 1),    # A^CTAGT
    "BglII": ("AGATCT", 1),   # A^GATCT
    "NcoI": ("CCATGG", 2),    # C^CATGG
    "NdeI": ("CATATG", 2),    # CA^TATG
    "XhoI": ("CTCGAG", 1),    # C^TCGAG
    "SalI": ("GTCGAC", 1),    # G^TCGAC
    "ClaI": ("ATCGAT", 2),    # AT^CGAT
    "HpaI": ("GTTAAC", 3),    # GTT^AAC
    "MluI": ("ACGCGT", 1),    # A^CGCGT
    "NheI": ("GCTAGC", 1),    # G^CTAGC
}

DNA_LADDER_100BP = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1200, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000
]


def parse_fasta(filepath: str) -> List[Tuple[str, str]]:
    """Parse a FASTA file and return list of (header, sequence)."""
    sequences = []
    current_header = None
    current_seq_parts = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header is not None:
                        sequences.append((current_header, ''.join(current_seq_parts)))
                    current_header = line[1:].strip()
                    current_seq_parts = []
                else:
                    if current_header is None:
                        raise ValueError("FASTA file missing header line")
                    if not re.match('^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$', line):
                        raise ValueError(f"Invalid DNA characters in line: {line[:50]}...")
                    current_seq_parts.append(line.upper())
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq_parts)))
        elif not sequences:
            raise ValueError("No sequences found in FASTA file")

        return sequences
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error parsing FASTA: {e}", file=sys.stderr)
        sys.exit(1)


def parse_enzyme_spec(spec: str) -> Tuple[str, str, int]:
    """Parse enzyme specification: NAME[:RECOG[:OFFSET]].
    If only name, look up in DB.
    Returns (name, recognition_seq, cut_offset).
    """
    if spec.upper() in ENZYME_DB:
        recog, offset = ENZYME_DB[spec.upper()]
        return (spec.upper(), recog, offset)

    parts = spec.split(':')
    if len(parts) == 1:
        raise ValueError(f"Unknown enzyme '{spec}'. Provide custom as NAME:SEQ:OFFSET.")
    if len(parts) == 2:
        recog = parts[1].upper()
        raise ValueError(f"Cut offset required for custom enzyme '{spec}'. Use NAME:SEQ:OFFSET.")
    if len(parts) == 3:
        name = parts[0].upper()
        recog = parts[1].upper()
        try:
            offset = int(parts[2])
        except ValueError:
            raise ValueError(f"Offset must be integer, got '{parts[2]}'")
        if not (0 <= offset <= len(recog)):
            raise ValueError(f"Offset {offset} out of range for recognition seq length {len(recog)}")
        return (name, recog, offset)
    raise ValueError(f"Malformed enzyme spec: '{spec}'")


def find_cut_positions(sequence: str, recognition: str, offset: int) -> List[int]:
    """Return all cut positions (0-based index after which the cut occurs) for a given enzyme."""
    cuts = []
    recog_len = len(recognition)
    start = 0
    while True:
        idx = sequence.find(recognition, start)
        if idx == -1:
            break
        cuts.append(idx + offset)
        start = idx + 1  # allow overlapping sites
    return sorted(set(cuts))


def digest_linear(seq_len: int, cut_positions: List[int], min_fragment: int) -> List[int]:
    """Compute fragment lengths for linear DNA."""
    all_cuts = sorted(cut_positions)
    fragments = []
    prev = 0
    for cut in all_cuts:
        frag_len = cut - prev
        if frag_len >= min_fragment:
            fragments.append(frag_len)
        prev = cut
    final_frag = seq_len - prev
    if final_frag >= min_fragment:
        fragments.append(final_frag)
    return fragments


def digest_circular(seq_len: int, cut_positions: List[int], min_fragment: int) -> List[int]:
    """Compute fragment lengths for circular DNA."""
    if not cut_positions:
        return [seq_len] if seq_len >= min_fragment else []
    sorted_cuts = sorted(cut_positions)
    fragments = []
    n = len(sorted_cuts)
    for i in range(n):
        start = sorted_cuts[i]
        end = sorted_cuts[(i + 1) % n]
        if i == n - 1:
            dist = seq_len - start + end
        else:
            dist = end - start
        if dist >= min_fragment:
            fragments.append(dist)
    return fragments


def print_fragment_table(fragments: List[int], enzyme_names: List[str]):
    """Print a simple table of fragments."""
    print("\n--- Fragment Summary ---")
    fragments_sorted = sorted(fragments, reverse=True)
    total_len = sum(fragments_sorted)
    print(f"{'#':>3} {'Size (bp)':>9}")
    print("-" * 14)
    for i, size in enumerate(fragments_sorted, 1):
        print(f"{i:>3} {size:>9}")
    print("-" * 14)
    print(f"{'Total':>3} {total_len:>9}")
    print()


def draw_ascii_gel(fragments: List[int], ladder_sizes: List[int], height: int = 30):
    """Draw ASCII gel with ladder and sample lane."""
    all_sizes = ladder_sizes + fragments
    if not all_sizes:
        print("No fragments to display.")
        return
    max_s = max(all_sizes)
    min_s = min(all_sizes)
    max_log = math.log10(max(max_s, 1))
    min_log = math.log10(max(min_s, 1))
    if max_log == min_log:
        max_log += 0.1  # avoid division by zero

    def row_for(size):
        log_s = math.log10(size) if size > 0 else min_log
        frac = (max_log - log_s) / (max_log - min_log)  # larger -> top (row 0)
        return max(0, min(height - 1, int(round(frac * (height - 1)))))

    # Build ladder row map: row -> list of ladder sizes
    ladder_rows = {}
    for sz in ladder_sizes:
        r = row_for(sz)
        ladder_rows.setdefault(r, []).append(sz)

    sample_bands = set()
    for sz in fragments:
        sample_bands.add(row_for(sz))

    print(f"{'Ladder (bp)':>12} | {'Sample':^5}")
    print("-" * 25)
    for r in range(height):
        labels = ",".join(str(s) for s in sorted(ladder_rows.get(r, []), reverse=True))
        band_char = 'X' if r in sample_bands else '|'
        line = f"{labels:>12} | {band_char:^5}"
        print(line)


def main():
    parser = argparse.ArgumentParser(
        description="Restriction Enzyme Digest Simulator - in silico restriction digestion of DNA."
    )
    parser.add_argument('--fasta', '-f', required=True, help='Path to FASTA file')
    parser.add_argument('--enzymes', '-e', required=True,
                        help='Comma-separated enzyme names or custom specs (NAME:SEQ:OFFSET)')
    parser.add_argument('--circular', '-c', action='store_true',
                        help='Treat DNA as circular (default: linear)')
    parser.add_argument('--min-fragment', '-m', type=int, default=1,
                        help='Minimum fragment length to report (default: 1)')
    parser.add_argument('--output', '-o', choices=['table', 'gel', 'both'], default='table',
                        help='Output format: table (default), gel, or both')
    parser.add_argument('--gel-height', type=int, default=30,
                        help='Height of ASCII gel (default: 30)')

    args = parser.parse_args()

    # Parse FASTA
    sequences = parse_fasta(args.fasta)

    # Parse enzymes
    enzyme_specs = [e.strip() for e in args.enzymes.split(',') if e.strip()]
    if not enzyme_specs:
        print("Error: No enzymes specified.", file=sys.stderr)
        sys.exit(1)

    enzymes = []
    enzyme_names = []
    for spec in enzyme_specs:
        try:
            name, recog, offset = parse_enzyme_spec(spec)
            enzymes.append((name, recog, offset))
            enzyme_names.append(name)
        except ValueError as e:
            print(f"Error parsing enzyme spec '{spec}': {e}", file=sys.stderr)
            sys.exit(1)

    # Process each sequence
    for header, seq in sequences:
        print(f"\n=== Digest of {header} ===")
        print(f"Sequence length: {len(seq)} bp, {'circular' if args.circular else 'linear'}")
        print(f"Enzymes: {', '.join(enzyme_names)}")

        all_cuts = []
        for name, recog, offset in enzymes:
            cuts = find_cut_positions(seq, recog, offset)
            if cuts:
                print(f"  {name} ({recog}) cut sites at: {cuts}")
            else:
                print(f"  {name} ({recog}): no cut sites found")
            all_cuts.extend(cuts)

        all_cuts = sorted(set(all_cuts))

        if args.circular:
            fragments = digest_circular(len(seq), all_cuts, args.min_fragment)
        else:
            fragments = digest_linear(len(seq), all_cuts, args.min_fragment)

        if args.output in ('table', 'both'):
            print_fragment_table(fragments, enzyme_names)
        if args.output in ('gel', 'both'):
            print("\n--- Simulated Gel Electrophoresis (100 bp ladder) ---")
            draw_ascii_gel(fragments, DNA_LADDER_100BP, height=args.gel_height)


if __name__ == '__main__':
    main()
