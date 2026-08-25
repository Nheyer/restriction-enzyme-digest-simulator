#!/usr/bin/env python3
"""Restriction Enzyme Digest Simulator."""

import argparse
import sys
import math
import re
from typing import Dict, List, NamedTuple, Tuple

from enzyme_data import (
    COMPLEMENT,
    DNA_LADDER_100BP,
    ENZYME_TABLE,
    IUPAC_ALPHABET,
    IUPAC_BASES,
)


class Enzyme(NamedTuple):
    """A restriction enzyme.

    ``cut_top`` and ``cut_bottom`` are 0-based offsets from the first base of
    the recognition site, both measured along the top strand: the cut falls
    immediately before the base at that offset. ``cut_bottom`` is where the
    enzyme nicks the complementary strand, expressed in the same top-strand
    frame, so that ``cut_bottom - cut_top`` is the overhang. Offsets may lie
    outside the site - type IIS enzymes such as BsaI cut downstream of it.
    """
    name: str
    site: str
    cut_top: int
    cut_bottom: int

    @property
    def is_palindromic(self) -> bool:
        """True if the site reads the same on both strands (EcoRI, not BsaI)."""
        return self.site == reverse_complement(self.site)

    @property
    def is_symmetric(self) -> bool:
        """True if binding in either orientation gives the same pair of cuts."""
        return self.is_palindromic and self.cut_bottom == len(self.site) - self.cut_top

    @property
    def overhang(self) -> int:
        """Positive for a 5' overhang, 0 for a blunt cut, negative for a 3' overhang."""
        return self.cut_bottom - self.cut_top

    @property
    def notation(self) -> str:
        """The enzyme's recognition specificity in NEB/REBASE notation."""
        return to_neb_notation(self)

    def describe_ends(self) -> str:
        if self.overhang > 0:
            return "5' overhang, {} nt".format(self.overhang)
        if self.overhang < 0:
            return "3' overhang, {} nt".format(-self.overhang)
        return "blunt"


# The static table lives in enzyme_data; this is it as typed Enzyme records.
ENZYME_DB: Dict[str, Enzyme] = {
    name: Enzyme(name, site, top, bottom)
    for name, (site, top, bottom) in ENZYME_TABLE.items()
}
# Enzyme names are matched case-insensitively but reported in canonical casing.
_ENZYME_LOOKUP = {name.upper(): name for name in ENZYME_DB}


def reverse_complement(seq: str) -> str:
    """Reverse complement of an IUPAC nucleotide string."""
    return seq.translate(COMPLEMENT)[::-1]


def validate_site(site: str) -> str:
    """Upper-case a recognition site and reject non-IUPAC characters."""
    site = site.strip().upper()
    if not site:
        raise ValueError("Recognition site is empty")
    bad = sorted(set(site) - set(IUPAC_ALPHABET))
    if bad:
        raise ValueError("Invalid IUPAC character(s) in recognition site: {}".format(', '.join(bad)))
    return site


def expand_ambiguity(site: str, max_expansions: int = 65536) -> List[str]:
    """Expand IUPAC ambiguity codes into every concrete ACGT sequence.

    HincII's GTYRAC becomes GTCAAC, GTCGAC, GTTAAC, GTTGAC. The count grows as
    the product of the codes' degeneracies, so ``max_expansions`` caps it:
    XcmI's nine Ns alone would give 262144 strings. Matching itself uses
    :func:`site_to_regex` and needs no expansion; this is for callers that want
    the explicit list.
    """
    site = validate_site(site)
    total = 1
    for code in site:
        total *= len(IUPAC_BASES[code])
        if total > max_expansions:
            raise ValueError(
                "Expanding '{}' would give more than {} sequences; raise "
                "max_expansions or match with site_to_regex() instead.".format(site, max_expansions)
            )
    results = ['']
    for code in site:
        results = [prefix + base for prefix in results for base in IUPAC_BASES[code]]
    return results


def site_to_regex(site: str) -> str:
    """Regex matching a recognition site, honouring IUPAC codes on both sides.

    A code in the site matches a base in the target only when that target base
    is unambiguously one of the bases the site accepts, so enzyme N matches a
    target N but enzyme A does not.
    """
    site = validate_site(site)
    parts = []
    for code in site:
        allowed = set(IUPAC_BASES[code])
        matching = [c for c in IUPAC_ALPHABET if set(IUPAC_BASES[c]) <= allowed]
        parts.append(matching[0] if len(matching) == 1 else '[' + ''.join(matching) + ']')
    return ''.join(parts)


def parse_neb_notation(text: str) -> Tuple[str, int, int]:
    """Parse a NEB/REBASE recognition specificity into (site, cut_top, cut_bottom).

    Accepts the two notations NEB publishes:
      * a caret inside the site - ``G^AATTC``, ``GGTAC^C``, ``^GATC``
      * cut offsets after the site - ``GGTCTC(1/5)``, ``GAATTC(-5/-1)``
    Enzymes that cut on both sides of their site (``(10/15)ACNNNNGTAYC(12/7)``)
    are rejected: two cuts per site do not fit the single-cut fragment model.
    """
    text = text.strip().upper()
    if not text:
        raise ValueError("Empty recognition specificity")
    if text.startswith('('):
        raise ValueError(
            "'{}' cuts on both sides of its recognition site; dual-cut enzymes "
            "(BaeI, BcgI, CspCI ...) are not supported.".format(text)
        )

    offsets = re.fullmatch(r'([^()]+)\((-?\d+)/(-?\d+)\)', text)
    if offsets:
        if '^' in offsets.group(1):
            raise ValueError("'{}' mixes caret and offset notation; use one or the other".format(text))
        site = validate_site(offsets.group(1))
        length = len(site)
        return site, length + int(offsets.group(2)), length + int(offsets.group(3))

    if '^' in text:
        if text.count('^') > 1:
            raise ValueError("'{}' has more than one cut position marker '^'".format(text))
        cut_top = text.index('^')
        site = validate_site(text.replace('^', ''))
        # A caret marks only the top-strand cut; the bottom cut is its mirror.
        return site, cut_top, len(site) - cut_top

    validate_site(text)  # surface bad characters before the generic message
    raise ValueError(
        "No cut position in '{}'. Mark it with '^' (G^AATTC), append offsets "
        "(GGTCTC(1/5)), or use NAME:SEQ:OFFSET.".format(text)
    )


def to_neb_notation(enzyme: Enzyme) -> str:
    """Render an enzyme as a NEB/REBASE recognition specificity string."""
    length = len(enzyme.site)
    if 0 <= enzyme.cut_top <= length and enzyme.cut_bottom == length - enzyme.cut_top:
        return enzyme.site[:enzyme.cut_top] + '^' + enzyme.site[enzyme.cut_top:]
    return "{}({}/{})".format(enzyme.site, enzyme.cut_top - length, enzyme.cut_bottom - length)


def to_legacy_spec(enzyme: Enzyme) -> str:
    """Render an enzyme as the NAME:SEQ:OFFSET spec this tool started with.

    That form records only the top-strand cut and assumes the bottom-strand cut
    is its mirror image, so it cannot represent an asymmetric cutter. Rather
    than silently dropping ``cut_bottom`` it raises for enzymes it would garble.
    """
    length = len(enzyme.site)
    if not 0 <= enzyme.cut_top <= length:
        raise ValueError(
            "{} cuts outside its recognition site ({}); NAME:SEQ:OFFSET cannot express "
            "that. Use NAME:{} instead.".format(enzyme.name, enzyme.notation, enzyme.notation)
        )
    if enzyme.cut_bottom != length - enzyme.cut_top:
        raise ValueError(
            "{} cuts the two strands asymmetrically ({}); NAME:SEQ:OFFSET would lose the "
            "bottom-strand cut. Use NAME:{} instead.".format(
                enzyme.name, enzyme.notation, enzyme.notation)
        )
    return "{}:{}:{}".format(enzyme.name, enzyme.site, enzyme.cut_top)


def parse_enzyme_spec(spec: str) -> Enzyme:
    """Parse an enzyme specification into an :class:`Enzyme`.

    Accepted forms:
      * ``EcoRI``               - a built-in, matched case-insensitively
      * ``NAME:G^AATTC``        - NEB notation, caret form
      * ``NAME:GGTCTC(1/5)``    - NEB notation, offset form (non-palindromic)
      * ``NAME:SEQ:OFFSET``     - top-strand cut, bottom strand mirrored
      * ``NAME:SEQ:TOP:BOTTOM`` - both cuts, as offsets from the site start
    """
    spec = spec.strip()
    canonical = _ENZYME_LOOKUP.get(spec.upper())
    if canonical:
        return ENZYME_DB[canonical]

    parts = spec.split(':')
    if len(parts) == 1:
        if '^' in spec or '(' in spec:
            # A bare recognition specificity, e.g. G^AATTC or GGTCTC(1/5).
            site, cut_top, cut_bottom = parse_neb_notation(spec)
            return Enzyme("custom", site, cut_top, cut_bottom)
        raise ValueError(
            "Unknown enzyme '{}'. Provide custom as NAME:SEQ:OFFSET or "
            "NAME:RECOGNITION_SITE (e.g. MyI:G^AATTC, MyI:GGTCTC(1/5)).".format(spec)
        )

    name = parts[0].strip()
    if not name:
        raise ValueError("Malformed enzyme spec: '{}' (missing enzyme name)".format(spec))

    if len(parts) == 2:
        recog = parts[1].strip()
        if '^' not in recog and '(' not in recog:
            raise ValueError(
                "Cut offset required for custom enzyme '{}'. Use NAME:SEQ:OFFSET or mark "
                "the cut in the site itself (NAME:{}^{}).".format(spec, recog[:1], recog[1:])
            )
        site, cut_top, cut_bottom = parse_neb_notation(recog)
        return Enzyme(name, site, cut_top, cut_bottom)

    if len(parts) in (3, 4):
        site = validate_site(parts[1])
        raw_offsets = parts[2:]
        try:
            cuts = [int(value) for value in raw_offsets]
        except ValueError:
            bad = next(v for v in raw_offsets if not _is_int(v))
            raise ValueError("Offset must be integer, got '{}'".format(bad))
        cut_top = cuts[0]
        cut_bottom = cuts[1] if len(cuts) == 2 else len(site) - cut_top
        if len(cuts) == 1 and not 0 <= cut_top <= len(site):
            raise ValueError(
                "Offset {} out of range for recognition seq length {}".format(cut_top, len(site))
            )
        return Enzyme(name, site, cut_top, cut_bottom)

    raise ValueError("Malformed enzyme spec: '{}'".format(spec))


def _is_int(value: str) -> bool:
    try:
        int(value)
        return True
    except ValueError:
        return False


class SiteHit(NamedTuple):
    """One recognition site found in a sequence.

    ``start`` is where the site begins on the top strand, ``strand`` is '+'
    when the site reads in the given sequence's direction and '-' when it was
    found as the reverse complement, and ``cut`` is the top-strand cut position
    in sequence coordinates, before normalisation to the sequence bounds.
    """
    start: int
    strand: str
    cut: int


def find_sites(sequence: str, enzyme: Enzyme) -> List[SiteHit]:
    """Find every recognition site for an enzyme, on both strands.

    A non-palindromic enzyme such as BsaI recognises its site on either strand,
    and a site on the bottom strand cuts the top strand on the other side of
    the site, so both orientations are searched. For a fully symmetric enzyme
    the two searches coincide, so each site is reported once.
    """
    length = len(enzyme.site)
    hits: Dict[Tuple[int, int], SiteHit] = {}
    orientations = [('+', enzyme.site)]
    if not enzyme.is_symmetric:
        orientations.append(('-', reverse_complement(enzyme.site)))

    for strand, pattern in orientations:
        # Lookahead, so that overlapping sites are all reported.
        regex = re.compile('(?=(' + site_to_regex(pattern) + '))')
        for match in regex.finditer(sequence):
            start = match.start()
            if strand == '+':
                cut = start + enzyme.cut_top
            else:
                # Mirror the enzyme frame: an offset p sits at start + length - p,
                # so the enzyme's bottom-strand cut lands on our top strand.
                cut = start + length - enzyme.cut_bottom
            hits.setdefault((start, cut), SiteHit(start, strand, cut))
    return sorted(hits.values())


def normalise_cuts(cuts: List[int], seq_len: int, circular: bool) -> Tuple[List[int], int]:
    """Bring cut positions into sequence coordinates.

    Type IIS enzymes cut a fixed distance away from their site, which can fall
    past the end - or before the start - of the sequence. On circular DNA that
    wraps around; on linear DNA there is nothing there to cut, so the cut is
    dropped. Returns the usable cut positions and how many were dropped.
    """
    kept, dropped = set(), 0
    for cut in cuts:
        if circular:
            if seq_len:
                kept.add(cut % seq_len)
        elif 1 <= cut <= seq_len - 1:
            kept.add(cut)
        else:
            dropped += 1
    return sorted(kept), dropped


def find_cut_positions(sequence: str, recognition: str, offset: int) -> List[int]:
    """Return top-strand cut positions for a recognition site and cut offset.

    Kept for callers using the original top-strand-only interface; it now
    understands IUPAC ambiguity codes. :func:`find_sites` also searches the
    bottom strand, which is what non-palindromic enzymes need.
    """
    regex = re.compile('(?=(' + site_to_regex(recognition) + '))')
    return sorted({match.start() + offset for match in regex.finditer(sequence)})


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


def print_enzyme_list():
    """Print the built-in enzymes with their recognition specificities."""
    print(f"{'Enzyme':<9} {'Recognition site':<22} {'Ends':<18} Type")
    print("-" * 64)
    for name in sorted(ENZYME_DB, key=str.lower):
        enzyme = ENZYME_DB[name]
        kind = "palindromic" if enzyme.is_palindromic else "non-palindromic"
        if not 0 <= enzyme.cut_top <= len(enzyme.site):
            kind += ", cuts outside site"
        print(f"{name:<9} {enzyme.notation:<22} {enzyme.describe_ends():<18} {kind}")
    print(f"\n{len(ENZYME_DB)} enzymes. Specificities from REBASE (rebase.neb.com).")


def print_conversion(spec: str):
    """Show an enzyme spec in both notations."""
    enzyme = parse_enzyme_spec(spec)
    print(f"Enzyme:            {enzyme.name}")
    print(f"Recognition site:  {enzyme.notation}")
    try:
        print(f"NAME:SEQ:OFFSET:   {to_legacy_spec(enzyme)}")
    except ValueError as e:
        print(f"NAME:SEQ:OFFSET:   not representable - {e}")
    print(f"Cut (top/bottom):  {enzyme.cut_top}/{enzyme.cut_bottom} from the start of the site")
    print(f"Ends:              {enzyme.describe_ends()}")
    print(f"Palindromic:       {'yes' if enzyme.is_palindromic else 'no'}")
    if set(enzyme.site) - set('ACGT'):
        try:
            variants = expand_ambiguity(enzyme.site)
            shown = ', '.join(variants[:8]) + (', ...' if len(variants) > 8 else '')
            print(f"Ambiguity codes:   {len(variants)} concrete sites ({shown})")
        except ValueError as e:
            print(f"Ambiguity codes:   too degenerate to expand - {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Restriction Enzyme Digest Simulator - in silico restriction digestion of DNA."
    )
    parser.add_argument('--fasta', '-f', help='Path to FASTA file')
    parser.add_argument('--enzymes', '-e',
                        help='Comma-separated enzyme names or custom specs '
                             '(NAME:SEQ:OFFSET, NAME:G^AATTC, NAME:GGTCTC(1/5))')
    parser.add_argument('--circular', '-c', action='store_true',
                        help='Treat DNA as circular (default: linear)')
    parser.add_argument('--min-fragment', '-m', type=int, default=1,
                        help='Minimum fragment length to report (default: 1)')
    parser.add_argument('--output', '-o', choices=['table', 'gel', 'both'], default='table',
                        help='Output format: table (default), gel, or both')
    parser.add_argument('--gel-height', type=int, default=30,
                        help='Height of ASCII gel (default: 30)')
    parser.add_argument('--list-enzymes', action='store_true',
                        help='List the built-in enzymes and exit')
    parser.add_argument('--convert', metavar='SPEC',
                        help='Show an enzyme in both notations and exit')

    args = parser.parse_args()

    if args.list_enzymes:
        print_enzyme_list()
        return

    if args.convert:
        try:
            print_conversion(args.convert)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        return

    if not args.fasta or not args.enzymes:
        parser.error("the following arguments are required: --fasta/-f, --enzymes/-e")

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
            enzyme = parse_enzyme_spec(spec)
            enzymes.append(enzyme)
            enzyme_names.append(enzyme.name)
        except ValueError as e:
            print(f"Error parsing enzyme spec '{spec}': {e}", file=sys.stderr)
            sys.exit(1)

    # Process each sequence
    for header, seq in sequences:
        print(f"\n=== Digest of {header} ===")
        print(f"Sequence length: {len(seq)} bp, {'circular' if args.circular else 'linear'}")
        print(f"Enzymes: {', '.join(enzyme_names)}")

        all_cuts = []
        for enzyme in enzymes:
            hits = find_sites(seq, enzyme)
            cuts, dropped = normalise_cuts([h.cut for h in hits], len(seq), args.circular)
            label = f"{enzyme.name} ({enzyme.notation})"
            if cuts:
                print(f"  {label} cut sites at: {cuts}")
            elif hits:
                print(f"  {label}: {len(hits)} recognition site(s) found, none cut")
            else:
                print(f"  {label}: no cut sites found")
            if hits and not enzyme.is_symmetric:
                bottom = sum(1 for h in hits if h.strand == '-')
                print(f"    {len(hits) - bottom} site(s) on the top strand, "
                      f"{bottom} on the bottom strand")
            if dropped:
                print(f"    {dropped} site(s) recognised but cutting beyond the end of the "
                      f"linear sequence - not cut")
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
