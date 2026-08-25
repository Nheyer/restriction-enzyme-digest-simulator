#!/usr/bin/env python3
"""Regenerate the enzyme tables in enzyme_data.py from restriction_enzymes.csv.

The CSV is the source of truth: a two-column ``enzyme,recognition_sequence``
export of NEB's commercially available specificities, written in the notation
NEB uses on the bench, where ``^`` marks the top-strand cut and ``_`` the
bottom-strand cut::

    EcoRI,G^AATTC              cuts both strands inside the site
    AatII,G_ACGT^C             3' overhang, the two cuts are not mirrored
    AcuI,CTGAAGNNNNNNNNNNNNNN_NN^   type IIS, cuts downstream of the site

enzyme_data.py has to stay pure data - no functions, so no CSV parsing there -
which is why the table is generated into it rather than read at import time.
Generating also keeps the tool a pair of flat modules: setuptools ships
py-modules but not data files alongside them, so a CSV read at runtime would
not survive ``pip install``. Pass a different CSV at runtime with --enzyme-db.

Usage:  python tools/build_enzyme_table.py [csv] [--check]
        --check exits non-zero if enzyme_data.py is out of date.
"""

import argparse
import os
import sys
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from enzyme_digest import load_enzyme_csv  # noqa: E402  (after sys.path setup)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CSV = os.path.join(ROOT, 'restriction_enzymes.csv')
TARGET = os.path.join(ROOT, 'enzyme_data.py')

BEGIN = "# --- BEGIN GENERATED TABLES (tools/build_enzyme_table.py) ---"
END = "# --- END GENERATED TABLES ---"


def render(standard, nicking, dual, csv_name) -> str:
    width = max(len(n) for n in standard) + 3
    lines = [BEGIN,
             "# Generated from {} - do not edit by hand.".format(csv_name),
             "# Regenerate with: python tools/build_enzyme_table.py",
             "#",
             "# (site, cut_top, cut_bottom): 0-based offsets from the first base of the",
             "# site, both measured along the top strand, the cut falling immediately",
             "# before the base at that offset. Offsets may sit outside the site - type",
             "# IIS enzymes cut downstream of it. The comment is the NEB specificity.",
             "ENZYME_TABLE: Dict[str, Tuple[str, int, int]] = {"]
    for name in sorted(standard, key=str.lower):
        site, top, bottom = standard[name]
        entry = '    {:<{w}} ({!r}, {}, {}),'.format('"%s":' % name, site, top, bottom, w=width)
        lines.append('{:<62} # {}'.format(entry, _notation(name, standard)))
    lines.append("}")
    lines.append("")
    lines.append("# Nicking enzymes cut one strand only. They leave the duplex intact, so")
    lines.append("# there is no fragment pattern to simulate; named here to give a better")
    lines.append("# error than 'unknown enzyme'.")
    lines.append("NICKING_ENZYMES: Dict[str, str] = {")
    for name in sorted(nicking, key=str.lower):
        lines.append('    "{}": "{}",'.format(name, nicking[name]))
    lines.append("}")
    lines.append("")
    lines.append("# Type IIB enzymes cut on both sides of their site, excising it. Two cuts")
    lines.append("# per site do not fit the single-cut Enzyme record, so they are listed")
    lines.append("# rather than loaded.")
    lines.append("DUAL_CUT_ENZYMES: Dict[str, str] = {")
    for name in sorted(dual, key=str.lower):
        lines.append('    "{}": "{}",'.format(name, dual[name]))
    lines.append("}")
    lines.append(END)
    return '\n'.join(lines)


def _notation(name, standard):
    site, top, bottom = standard[name]
    length = len(site)
    if 0 <= top <= length and bottom == length - top:
        return site[:top] + '^' + site[top:]
    return "{}({}/{})".format(site, top - length, bottom - length)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('csv', nargs='?', default=DEFAULT_CSV)
    parser.add_argument('--check', action='store_true',
                        help='exit non-zero if enzyme_data.py is out of date')
    args = parser.parse_args()

    tables = load_enzyme_csv(args.csv)
    standard = {name: (e.site, e.cut_top, e.cut_bottom)
                for name, e in tables.enzymes.items()}
    nicking, dual = tables.nicking, tables.dual_cut

    block = render(standard, nicking, dual, os.path.basename(args.csv))
    with open(TARGET, encoding='utf-8') as handle:
        current = handle.read()

    if BEGIN not in current or END not in current:
        raise SystemExit("markers not found in {}".format(TARGET))
    head, rest = current.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    updated = head + block + tail

    if args.check:
        if updated != current:
            print("enzyme_data.py is out of date; run tools/build_enzyme_table.py",
                  file=sys.stderr)
            return 1
        print("enzyme_data.py is up to date ({} enzymes)".format(len(standard)))
        return 0

    with open(TARGET, 'w', encoding='utf-8') as handle:
        handle.write(updated)
    print("wrote {}: {} standard, {} nicking, {} dual-cut".format(
        TARGET, len(standard), len(nicking), len(dual)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
