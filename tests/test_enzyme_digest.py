import ast
import importlib
import pathlib
import random
import re

import pytest
import tempfile
import os
import enzyme_data
import enzyme_digest
from enzyme_digest import (
    Enzyme,
    parse_fasta,
    parse_enzyme_spec,
    parse_neb_notation,
    to_neb_notation,
    to_legacy_spec,
    expand_ambiguity,
    site_to_regex,
    reverse_complement,
    find_sites,
    find_cut_positions,
    normalise_cuts,
    digest_linear,
    digest_circular,
    strip_gaps,
    sequence_alphabet,
    ENZYME_DB,
    AMBIGUITY_DEFINITE,
    AMBIGUITY_POSSIBLE,
)
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      'enzyme_digest.py')


def write_fasta(content):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(content)
        return f.name


def run_cli_with_status(content, *args):
    """Run the CLI over a temporary FASTA; returns (combined output, exit code)."""
    fname = write_fasta(content)
    try:
        proc = subprocess.run([sys.executable, SCRIPT, '-f', fname] + list(args),
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return proc.stdout.decode(), proc.returncode
    finally:
        os.unlink(fname)


def run_cli(content, *args):
    return run_cli_with_status(content, *args)[0]


class TestParseFasta:
    def test_basic(self):
        content = ">seq1\nATGC\nTAGG\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(content)
            fname = f.name
        try:
            seqs = parse_fasta(fname)
            assert len(seqs) == 1
            assert seqs[0][0] == "seq1"
            assert seqs[0][1] == "ATGCTAGG"
        finally:
            os.unlink(fname)

    def test_multiple_sequences(self):
        content = ">h1\nAAAA\n>h2\nGGGG\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(content)
            fname = f.name
        try:
            seqs = parse_fasta(fname)
            assert len(seqs) == 2
            assert seqs[1][0] == "h2"
            assert seqs[1][1] == "GGGG"
        finally:
            os.unlink(fname)

    def test_invalid_characters_exits(self):
        content = ">test\nATGC!XYZ\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(content)
            fname = f.name
        try:
            with pytest.raises(SystemExit):
                parse_fasta(fname)
        finally:
            os.unlink(fname)

    def test_no_header_exits(self):
        content = "ATGC\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(content)
            fname = f.name
        try:
            with pytest.raises(SystemExit):
                parse_fasta(fname)
        finally:
            os.unlink(fname)


class TestParseEnzymeSpec:
    def test_builtin(self):
        enzyme = parse_enzyme_spec("EcoRI")
        assert enzyme.name == "EcoRI"
        assert enzyme.site == "GAATTC"
        assert enzyme.cut_top == 1
        assert enzyme.cut_bottom == 5

    def test_case_insensitive_builtin(self):
        # Names are matched case-insensitively but reported canonically.
        assert parse_enzyme_spec("ecori").name == "EcoRI"
        assert parse_enzyme_spec("BAMHI") == ENZYME_DB["BamHI"]

    def test_custom_valid(self):
        enzyme = parse_enzyme_spec("MyI:ACG:2")
        assert enzyme.name == "MyI"
        assert enzyme.site == "ACG"
        assert enzyme.cut_top == 2
        # The legacy form mirrors the cut onto the bottom strand.
        assert enzyme.cut_bottom == 1

    def test_custom_missing_offset_raises(self):
        with pytest.raises(ValueError, match="offset required"):
            parse_enzyme_spec("Test:GAATTC")

    def test_custom_offset_out_of_range_raises(self):
        with pytest.raises(ValueError, match="out of range"):
            parse_enzyme_spec("Bad:ACG:5")

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError, match="Unknown enzyme"):
            parse_enzyme_spec("NonExistent")

    def test_bare_recognition_string_without_name(self):
        enzyme = parse_enzyme_spec("G^AATTC")
        assert (enzyme.site, enzyme.cut_top, enzyme.cut_bottom) == ("GAATTC", 1, 5)
        assert parse_enzyme_spec("GGTCTC(1/5)").cut_top == 7

    def test_bare_dual_cut_string_rejected(self):
        with pytest.raises(ValueError, match="dual-cut"):
            parse_enzyme_spec("(10/15)ACNNNNGTAYC(12/7)")

    def test_caret_notation(self):
        enzyme = parse_enzyme_spec("MyI:G^AATTC")
        assert (enzyme.site, enzyme.cut_top, enzyme.cut_bottom) == ("GAATTC", 1, 5)

    def test_offset_notation_non_palindromic(self):
        enzyme = parse_enzyme_spec("MyI:GGTCTC(1/5)")
        assert (enzyme.site, enzyme.cut_top, enzyme.cut_bottom) == ("GGTCTC", 7, 11)
        assert not enzyme.is_palindromic

    def test_explicit_two_strand_offsets(self):
        enzyme = parse_enzyme_spec("MyI:GGTCTC:7:11")
        assert enzyme == ENZYME_DB["BsaI"]._replace(name="MyI")

    def test_lowercase_site_is_normalised(self):
        assert parse_enzyme_spec("MyI:gaattc:1").site == "GAATTC"

    def test_invalid_site_character_raises(self):
        with pytest.raises(ValueError, match="Invalid IUPAC"):
            parse_enzyme_spec("MyI:GAXTTC:1")

    def test_non_integer_offset_raises(self):
        with pytest.raises(ValueError, match="Offset must be integer"):
            parse_enzyme_spec("MyI:GAATTC:x")


class TestNebNotation:
    @pytest.mark.parametrize("text,expected", [
        ("G^AATTC", ("GAATTC", 1, 5)),      # EcoRI, 5' overhang
        ("GGTAC^C", ("GGTACC", 5, 1)),      # KpnI, 3' overhang
        ("CCC^GGG", ("CCCGGG", 3, 3)),      # SmaI, blunt
        ("^GATC", ("GATC", 0, 4)),          # MboI, cuts before the site
        ("GGTCTC(1/5)", ("GGTCTC", 7, 11)),  # BsaI, type IIS
        ("GCTCTTC(1/4)", ("GCTCTTC", 8, 11)),  # SapI
        ("GAATTC(-5/-1)", ("GAATTC", 1, 5)),  # EcoRI in offset form
    ])
    def test_parse(self, text, expected):
        assert parse_neb_notation(text) == expected

    def test_round_trip_over_whole_db(self):
        for name, enzyme in ENZYME_DB.items():
            site, top, bottom = parse_neb_notation(to_neb_notation(enzyme))
            assert (site, top, bottom) == (enzyme.site, enzyme.cut_top, enzyme.cut_bottom), name

    def test_round_trip_through_spec_string(self):
        for name, enzyme in ENZYME_DB.items():
            assert parse_enzyme_spec(f"{name}:{enzyme.notation}") == enzyme, name

    def test_dual_cut_enzyme_rejected(self):
        with pytest.raises(ValueError, match="dual-cut"):
            parse_neb_notation("(10/15)ACNNNNGTAYC(12/7)")

    def test_two_carets_rejected(self):
        with pytest.raises(ValueError, match="more than one cut position"):
            parse_neb_notation("G^AAT^TC")

    def test_site_without_cut_rejected(self):
        with pytest.raises(ValueError, match="No cut position"):
            parse_neb_notation("GAATTC")


class TestLegacySpecConversion:
    def test_symmetric_enzyme_round_trips(self):
        for enzyme in ENZYME_DB.values():
            if not enzyme.is_symmetric:
                continue
            assert parse_enzyme_spec(to_legacy_spec(enzyme)) == enzyme, enzyme.name

    def test_asymmetric_enzyme_refuses_to_downgrade(self):
        # Every built-in that cuts inside its site mirrors the cut, so build one
        # that does not: NAME:SEQ:OFFSET would silently lose the bottom cut.
        with pytest.raises(ValueError, match="asymmetric"):
            to_legacy_spec(Enzyme("MyI", "GAATTC", 1, 3))

    def test_cut_outside_site_refuses_to_downgrade(self):
        with pytest.raises(ValueError, match="outside its recognition site"):
            to_legacy_spec(ENZYME_DB["BsaI"])


class TestPackaging:
    """The project has to stay installable, and the metadata has to match."""

    def test_version_is_a_release_number(self):
        assert re.fullmatch(r'\d+\.\d+(\.\d+)?([ab]\d+|rc\d+|\.dev\d+)?', enzyme_digest.__version__)

    def test_console_script_target_exists(self):
        # pyproject declares enzyme-digest = "enzyme_digest:main".
        assert callable(enzyme_digest.main)

    def test_pyproject_ships_every_module_the_tool_needs(self):
        """An installed copy has to carry enzyme_data too, not just the logic."""
        tomllib = pytest.importorskip("tomllib")  # 3.11+; skipped on older
        root = pathlib.Path(enzyme_digest.__file__).parent
        pyproject = root / "pyproject.toml"
        if not pyproject.exists():
            pytest.skip("running against an installed copy, not a checkout")
        config = tomllib.loads(pyproject.read_text(encoding='utf-8'))

        declared = set(config['tool']['setuptools']['py-modules'])
        not_shipped = {'setup', 'conftest'}  # tooling, not part of the package
        on_disk = {path.stem for path in root.glob('*.py')} - not_shipped
        assert on_disk <= declared, "module in the repo root missing from py-modules"

        entry = config['project']['scripts']['enzyme-digest']
        module, _, attribute = entry.partition(':')
        assert module in declared
        assert hasattr(importlib.import_module(module), attribute)

        version_attr = config['tool']['setuptools']['dynamic']['version']['attr']
        assert version_attr == 'enzyme_digest.__version__'


class TestDataModule:
    """The lookup tables live in enzyme_data, apart from the logic."""

    def test_table_and_db_agree(self):
        assert set(enzyme_data.ENZYME_TABLE) == set(ENZYME_DB)
        for name, (site, top, bottom) in enzyme_data.ENZYME_TABLE.items():
            assert ENZYME_DB[name] == Enzyme(name, site, top, bottom)

    def test_data_module_holds_no_logic(self):
        # No functions or classes of its own, so the tables cannot drift back
        # into logic.
        defined_here = [
            name for name, value in vars(enzyme_data).items()
            if getattr(value, '__module__', None) == 'enzyme_data'
        ]
        assert defined_here == []

    def test_data_module_does_not_import_the_logic(self):
        # It has to stay importable on its own; enzyme_digest depends on it,
        # never the other way round.
        tree = ast.parse(pathlib.Path(enzyme_data.__file__).read_text(encoding='utf-8'))
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.add(node.module)
        assert 'enzyme_digest' not in imported
        assert imported == {'typing'}

    def test_complement_table_is_self_inverse(self):
        for code in enzyme_data.IUPAC_ALPHABET:
            assert reverse_complement(reverse_complement(code)) == code
            # Complementing a code complements the bases it stands for.
            complemented = reverse_complement(code)
            assert set(enzyme_data.IUPAC_BASES[complemented]) == {
                {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[base]
                for base in enzyme_data.IUPAC_BASES[code]
            }


class TestEnzymeDatabase:
    def test_known_overhangs(self):
        # A typo in either offset shows up as a wrong overhang.
        assert ENZYME_DB["EcoRI"].overhang == 4      # 5' overhang, AATT
        assert ENZYME_DB["SmaI"].overhang == 0       # blunt
        assert ENZYME_DB["SacI"].overhang == -4      # 3' overhang
        assert ENZYME_DB["BsaI"].overhang == 4       # 5' overhang, 4 nt
        assert ENZYME_DB["MlyI"].overhang == 0       # blunt type IIS

    def test_overhang_is_self_complementary(self):
        """A palindromic enzyme cutting inside its site must leave a self-
        complementary overhang - the symmetry of the site forces it, so this
        catches an offset that is wrong in either strand.
        """
        checked = 0
        for name, enzyme in ENZYME_DB.items():
            inside = 0 <= enzyme.cut_top <= len(enzyme.site) and 0 <= enzyme.cut_bottom <= len(enzyme.site)
            if not enzyme.is_palindromic or not inside:
                continue
            checked += 1
            lo, hi = sorted((enzyme.cut_top, enzyme.cut_bottom))
            overhang = enzyme.site[lo:hi]
            assert overhang == reverse_complement(overhang), (name, enzyme.notation, overhang)
        assert checked > 60  # the check must actually be reaching the database

    def test_overhang_check_rejects_a_wrong_offset(self):
        # EcoRI with the bottom cut one base off leaves AAT/ATT, not AATT.
        wrong = Enzyme("EcoRI-bad", "GAATTC", 1, 4)
        overhang = wrong.site[1:4]
        assert overhang != reverse_complement(overhang)

    def test_sites_are_valid_iupac(self):
        for name, enzyme in ENZYME_DB.items():
            assert set(enzyme.site) <= set("ACGTRYSWKMBDHVN"), name
            assert enzyme.name == name

    def test_palindrome_classification(self):
        assert ENZYME_DB["EcoRI"].is_palindromic
        assert ENZYME_DB["BstXI"].is_palindromic     # CCANNNNNNTGG reads both ways
        assert ENZYME_DB["BstXI"].is_symmetric       # and its cuts mirror
        assert not ENZYME_DB["BsaI"].is_palindromic  # GGTCTC does not
        # Palindromic but asymmetric: both orientations have to be searched.
        assert Enzyme("MyI", "GAATTC", 1, 3).is_palindromic
        assert not Enzyme("MyI", "GAATTC", 1, 3).is_symmetric


class TestAmbiguityCodes:
    def test_expand_hincii(self):
        assert expand_ambiguity("GTYRAC") == ["GTCAAC", "GTCGAC", "GTTAAC", "GTTGAC"]

    def test_expand_unambiguous_site_is_itself(self):
        assert expand_ambiguity("GAATTC") == ["GAATTC"]

    def test_expansion_cap(self):
        with pytest.raises(ValueError, match="max_expansions"):
            expand_ambiguity("NNNNNNNNN", max_expansions=1000)

    def test_expansion_matches_regex_engine(self):
        # Every expanded string must be found by the regex built from the site.
        import re
        for site in ("GTYRAC", "CCWWGG", "GANTC", "RGCGCY"):
            regex = re.compile(site_to_regex(site))
            for variant in expand_ambiguity(site):
                assert regex.fullmatch(variant), (site, variant)

    def test_ambiguous_enzyme_finds_all_variants(self):
        seq = "AAA" + "GTCAAC" + "TTT" + "GTTGAC" + "AAA"
        cuts = [h.cut for h in find_sites(seq, ENZYME_DB["HincII"])]
        assert cuts == [6, 15]  # GTY^RAC in both variants

    def test_target_ambiguity_only_matches_when_unambiguous(self):
        # Enzyme N accepts a target N; enzyme A does not.
        assert find_sites("AAANNNNTTC", ENZYME_DB["XmnI"]) == []
        assert find_sites("GAANNNNTTC", ENZYME_DB["XmnI"]) != []
        assert find_cut_positions("GANTC", "GANTC", 1) == [1]
        assert find_cut_positions("GNNTC", "GANTC", 1) == []


class TestFindSites:
    def test_two_sites(self):
        seq = "GAATTCGAATTC"
        pos = find_cut_positions(seq, "GAATTC", offset=1)
        assert pos == [1, 7]

    def test_no_site(self):
        seq = "AAAA"
        pos = find_cut_positions(seq, "GAATTC", offset=1)
        assert pos == []

    def test_overlapping_sites(self):
        # "ATATAT" with enzyme recognition "ATAT" offset 2
        # positions 0-> cut at 2, then start 1 find again at 1, cut at 3
        seq = "ATATAT"
        pos = find_cut_positions(seq, "ATAT", offset=2)
        assert pos == [2, 4]

    def test_palindromic_site_reported_once(self):
        hits = find_sites("AAAGAATTCAAA", ENZYME_DB["EcoRI"])
        assert [(h.start, h.strand, h.cut) for h in hits] == [(3, '+', 4)]

    def test_non_palindromic_site_on_bottom_strand(self):
        bsai = ENZYME_DB["BsaI"]
        # GAGACC is the reverse complement of BsaI's GGTCTC.
        seq = "AAAAAAAAAA" + "GAGACC" + "AAAAAAAAAA"
        hits = find_sites(seq, bsai)
        assert [(h.start, h.strand, h.cut) for h in hits] == [(10, '-', 5)]

    def test_non_palindromic_site_on_both_strands(self):
        bsai = ENZYME_DB["BsaI"]
        seq = "AAAA" + "GGTCTC" + "A" * 20 + "GAGACC" + "AAAA"
        hits = find_sites(seq, bsai)
        assert sorted(h.strand for h in hits) == ['+', '-']
        top = next(h for h in hits if h.strand == '+')
        bottom = next(h for h in hits if h.strand == '-')
        assert top.cut == top.start + 7                    # 1 nt spacer downstream
        assert bottom.cut == bottom.start + 6 - 11         # mirrored, upstream

    def test_bottom_strand_hit_mirrors_top_strand_hit(self):
        # The same site read from the other strand must give the same geometry.
        bsai = ENZYME_DB["BsaI"]
        seq = "T" * 20 + "GGTCTC" + "T" * 20
        rc = reverse_complement(seq)
        forward = find_sites(seq, bsai)[0]
        reverse = find_sites(rc, bsai)[0]
        assert forward.strand == '+' and reverse.strand == '-'
        # Reverse-complementing swaps the strands: the cut the enzyme makes on
        # the bottom strand of `seq` is the one it makes on the top strand of `rc`.
        assert reverse.cut == len(seq) - (forward.start + bsai.cut_bottom)
        assert forward.cut == len(seq) - (reverse.start + len(bsai.site) - bsai.cut_top)


class TestNormaliseCuts:
    def test_linear_drops_cuts_outside_sequence(self):
        cuts, dropped = normalise_cuts([-3, 0, 5, 100, 101], 100, circular=False)
        assert cuts == [5]
        assert dropped == 4

    def test_circular_wraps_cuts(self):
        cuts, dropped = normalise_cuts([-2, 5, 103], 100, circular=True)
        assert cuts == [3, 5, 98]
        assert dropped == 0

    def test_type_iis_site_near_linear_end_is_not_cut(self):
        # BsaI cuts 7 nt past the site start, one base beyond this sequence.
        seq = "A" * 20 + "GGTCTC"
        hits = find_sites(seq, ENZYME_DB["BsaI"])
        assert len(hits) == 1
        cuts, dropped = normalise_cuts([h.cut for h in hits], len(seq), circular=False)
        assert cuts == []
        assert dropped == 1
        assert digest_linear(len(seq), cuts, 1) == [len(seq)]

    def test_type_iis_site_near_circular_end_wraps(self):
        seq = "A" * 20 + "GGTCTC" + "AAA"
        hits = find_sites(seq, ENZYME_DB["BsaI"])
        cuts, dropped = normalise_cuts([h.cut for h in hits], len(seq), circular=True)
        assert cuts == [(20 + 7) % len(seq)]
        assert dropped == 0

    def test_bottom_strand_site_at_sequence_start_is_not_cut(self):
        # A bottom-strand BsaI site here would cut 5 nt before the sequence.
        seq = "GAGACC" + "A" * 20
        hits = find_sites(seq, ENZYME_DB["BsaI"])
        assert [h.cut for h in hits] == [-5]
        cuts, dropped = normalise_cuts([h.cut for h in hits], len(seq), circular=False)
        assert (cuts, dropped) == ([], 1)
        assert all(f > 0 for f in digest_linear(len(seq), cuts, 1))


class TestDigestLinear:
    def test_basic(self):
        frags = digest_linear(10, [2, 5, 8], min_fragment=1)
        assert frags == [2, 3, 3, 2]

    def test_no_cuts(self):
        assert digest_linear(100, [], 1) == [100]

    def test_min_fragment_filter(self):
        frags = digest_linear(10, [2, 5, 8], min_fragment=3)
        assert frags == [3, 3]

    def test_end_fragments(self):
        # cuts at 0 not allowed (positions are >=1), but test
        frags = digest_linear(30, [10, 20], min_fragment=1)
        assert frags == [10, 10, 10]


class TestDigestCircular:
    def test_basic(self):
        frags = digest_circular(10, [2, 5, 8], min_fragment=1)
        assert sorted(frags) == [3, 3, 4]

    def test_no_cuts(self):
        assert digest_circular(100, [], 1) == [100]

    def test_single_cut(self):
        frags = digest_circular(50, [20], min_fragment=1)
        assert frags == [50]

    def test_min_fragment(self):
        frags = digest_circular(10, [2, 5, 8], min_fragment=4)
        assert frags == [4]


def _digest(seq, enzyme, circular):
    hits = find_sites(seq, enzyme)
    cuts, _ = normalise_cuts([h.cut for h in hits], len(seq), circular)
    if circular:
        return digest_circular(len(seq), cuts, 1)
    return digest_linear(len(seq), cuts, 1)


class TestFragmentInvariants:
    """Whatever the enzyme, the fragments must add back up to the sequence."""

    @pytest.mark.parametrize("circular", [False, True])
    def test_fragments_sum_to_sequence_length(self, circular):
        rng = random.Random(20260824)
        sequences = [''.join(rng.choice("ACGT") for _ in range(400)) for _ in range(3)]
        # Plus sequences seeded with sites right at both boundaries.
        for enzyme in ENZYME_DB.values():
            sequences.append(enzyme.site + ''.join(rng.choice("ACGT") for _ in range(30)))
            sequences.append(''.join(rng.choice("ACGT") for _ in range(30)) + enzyme.site)
            sequences.append(reverse_complement(enzyme.site) + "ACGT" * 5)
        for enzyme in ENZYME_DB.values():
            for seq in sequences:
                fragments = _digest(seq, enzyme, circular)
                assert sum(fragments) == len(seq), (enzyme.name, seq)
                assert all(f > 0 for f in fragments), (enzyme.name, seq)

    def test_every_enzyme_cuts_its_own_site(self):
        rng = random.Random(7)
        for enzyme in ENZYME_DB.values():
            flank = ''.join(rng.choice("ACGT") for _ in range(40))
            seq = flank + enzyme.site + flank
            hits = find_sites(seq, enzyme)
            assert any(h.start == 40 for h in hits), enzyme.name

    def test_enzyme_finds_site_on_reverse_complement_strand(self):
        rng = random.Random(11)
        for enzyme in ENZYME_DB.values():
            flank = ''.join(rng.choice("ACGT") for _ in range(40))
            seq = reverse_complement(flank + enzyme.site + flank)
            assert find_sites(seq, enzyme), enzyme.name


class TestEnzymeType:
    def test_notation_property_matches_helper(self):
        enzyme = Enzyme("MyI", "GAATTC", 1, 5)
        assert enzyme.notation == "G^AATTC" == to_neb_notation(enzyme)

    def test_describe_ends(self):
        assert Enzyme("a", "GAATTC", 1, 5).describe_ends() == "5' overhang, 4 nt"
        assert Enzyme("b", "CCCGGG", 3, 3).describe_ends() == "blunt"
        assert Enzyme("c", "GAGCTC", 5, 1).describe_ends() == "3' overhang, 4 nt"


class TestSequenceAmbiguity:
    """Ambiguity codes in the *target sequence* (the site side was already
    handled). 'definite' answers "where is this certain to cut?"; 'possible'
    answers "could this cut at all?" - the question behind checking that an
    enzyme leaves a sequence intact."""

    SSPI = ENZYME_DB["SspI"]

    def _cuts(self, seq, enzyme, mode, circular=False):
        return sorted({h.cut for h in find_sites(seq, enzyme, mode, circular)})

    def test_ambiguous_target_is_not_a_definite_cut(self):
        # SspI AAT^ATT vs GGAAYATYATR: the Ys could be C, so no guarantee.
        assert self._cuts("GGAAYATYATR", self.SSPI, AMBIGUITY_DEFINITE) == []

    def test_ambiguous_target_is_a_possible_cut(self):
        assert self._cuts("GGAAYATYATR", self.SSPI, AMBIGUITY_POSSIBLE) == [5]

    def test_n_in_target_is_possible_but_not_definite(self):
        eco = ENZYME_DB["EcoRI"]
        assert self._cuts("NAATTC", eco, AMBIGUITY_DEFINITE) == []
        assert self._cuts("NAATTC", eco, AMBIGUITY_POSSIBLE) == [1]

    def test_disjoint_codes_never_match(self):
        # Enzyme A against target Y (C or T) shares no base, in either mode.
        for mode in (AMBIGUITY_DEFINITE, AMBIGUITY_POSSIBLE):
            assert re.fullmatch(site_to_regex('A', mode), 'Y') is None

    def test_overlapping_codes_match_only_in_possible_mode(self):
        # Enzyme A against target R (A or G) may resolve to the A.
        assert re.fullmatch(site_to_regex('A', AMBIGUITY_DEFINITE), 'R') is None
        assert re.fullmatch(site_to_regex('A', AMBIGUITY_POSSIBLE), 'R')

    def test_possible_is_always_a_superset_of_definite(self):
        rng = random.Random(17)
        names = list(ENZYME_DB)
        for _ in range(400):
            seq = ''.join(rng.choice("ACGTRYKMSWBDHVN") for _ in range(rng.randint(8, 40)))
            enzyme = ENZYME_DB[rng.choice(names)]
            circ = rng.choice([True, False])
            definite = set(self._cuts(seq, enzyme, AMBIGUITY_DEFINITE, circ))
            possible = set(self._cuts(seq, enzyme, AMBIGUITY_POSSIBLE, circ))
            assert definite <= possible, (seq, enzyme.name)

    def test_unknown_mode_raises(self):
        with pytest.raises(ValueError, match="Unknown ambiguity mode"):
            site_to_regex("GAATTC", "maybe")


class TestCircularOriginSpanning:
    def test_site_spanning_the_origin_is_found(self):
        # GAATTC straddles the join: "TTC" at the end + "GAA" at the start.
        seq = "TTCGGGGGGGGGGGAA"
        eco = ENZYME_DB["EcoRI"]
        assert [h.cut for h in find_sites(seq, eco)] == []
        cuts, _ = normalise_cuts(
            [h.cut for h in find_sites(seq, eco, AMBIGUITY_DEFINITE, True)], len(seq), True)
        assert cuts == [14]

    def test_no_spurious_sites_when_nothing_spans_the_origin(self):
        seq = "AAGAATTCTTTTTTTTTT"
        eco = ENZYME_DB["EcoRI"]
        linear = [h.cut for h in find_sites(seq, eco)]
        circular = [h.cut for h in find_sites(seq, eco, AMBIGUITY_DEFINITE, True)]
        assert linear == circular

    def test_tandem_sites_each_counted_once(self):
        seq = "GAATTC" * 3
        cuts, _ = normalise_cuts(
            [h.cut for h in find_sites(seq, ENZYME_DB["EcoRI"], AMBIGUITY_DEFINITE, True)],
            len(seq), True)
        assert cuts == [1, 7, 13]

    def test_site_longer_than_sequence(self):
        assert find_sites("AAG", ENZYME_DB["EcoRI"], AMBIGUITY_DEFINITE, True) == []


class TestSequenceAlphabet:
    def test_strict_alphabet_rejects_ambiguity_codes(self):
        fname = write_fasta(">s\nACGTY\n")
        try:
            with pytest.raises(SystemExit):
                parse_fasta(fname, allow_ambiguity=False)
        finally:
            os.unlink(fname)

    def test_permissive_alphabet_accepts_iupac_and_gaps(self):
        fname = write_fasta(">s\nACGTRYKMSWBDHVN-\n")
        try:
            assert parse_fasta(fname, allow_ambiguity=True)[0][1] == "ACGTRYKMSWBDHVN-"
        finally:
            os.unlink(fname)

    def test_permissive_alphabet_still_rejects_junk(self):
        fname = write_fasta(">s\nACGTU\n")
        try:
            with pytest.raises(SystemExit):
                parse_fasta(fname, allow_ambiguity=True)
        finally:
            os.unlink(fname)

    def test_alphabet_contents(self):
        assert sequence_alphabet(False) == frozenset("ACGT")
        assert sequence_alphabet(True) == frozenset("ACGTRYKMSWBDHVN-")


class TestStripGaps:
    def test_removes_gaps_and_counts_them(self):
        assert strip_gaps("AAT-ATT") == ("AATATT", 1)

    def test_no_gaps_is_a_passthrough(self):
        assert strip_gaps("AATATT") == ("AATATT", 0)

    def test_gapped_site_is_a_real_site(self):
        # A gap is an alignment artefact; the molecule reads straight through.
        seq, removed = strip_gaps("GG-AAT--ATTGG")
        assert removed == 3
        assert [h.cut for h in find_sites(seq, ENZYME_DB["SspI"])] == [5]


class TestRejectionMessages:
    """Rejecting bad input only helps if the message says what is wrong and what
    to do about it, so assert on content, not just exit status."""

    def test_ambiguity_without_flag_is_rejected_loudly(self):
        out, code = run_cli_with_status(">plasmid_v2\nGGAAYATYATRNN\n", '-e', 'SspI')
        assert code != 0
        for char in ("'N'", "'R'", "'Y'"):
            assert char in out
        assert "plasmid_v2" in out and "line 2" in out
        assert "--ambiguity definite" in out and "--ambiguity possible" in out

    def test_junk_with_flag_is_rejected_and_alphabet_listed(self):
        out, code = run_cli_with_status(">bad\nACGTUZ\n", '-e', 'SspI', '-a', 'possible')
        assert code != 0
        assert "'U'" in out and "'Z'" in out and "ACGTRYKMSWBDHVN" in out
        assert "rerun with" not in out

    def test_mixed_junk_and_codes_are_distinguished(self):
        out, code = run_cli_with_status(">mixed\nACGTYU\n", '-e', 'SspI')
        assert "'U' is not valid DNA under any setting" in out
        assert "'Y' is an IUPAC ambiguity/gap code" in out

    def test_offending_character_is_never_truncated_away(self):
        out, code = run_cli_with_status(">long\n" + "A" * 200 + "Y\n", '-e', 'SspI')
        assert code != 0 and "'Y'" in out

    def test_valid_input_still_succeeds(self):
        out, code = run_cli_with_status(">ok\nAAGAATTCTT\n", '-e', 'EcoRI')
        assert code == 0 and "REJECTED" not in out


class TestAmbiguityCLI:
    def test_definite_mode_hints_at_hidden_sites(self):
        out = run_cli(">s\nGGAAYATYATR\n", '-e', 'SspI', '-a', 'definite')
        assert "no cut sites found" in out
        assert "1 further site(s) could cut" in out

    def test_possible_mode_marks_speculative_cuts(self):
        out = run_cli(">s\nGGAAYATYATR\n", '-e', 'SspI', '-a', 'possible')
        assert "cut sites at: [5?]" in out
        assert "maximal-cut scenario" in out

    def test_linear_run_flags_that_it_did_not_wrap(self):
        out = run_cli(">s\nGGAAYATYATR\n", '-e', 'BamHI', '-a', 'possible')
        assert "no cut possible under any resolution" in out
        assert "pass --circular" in out

    def test_circular_run_makes_the_claim_unqualified(self):
        out = run_cli(">s\nGGAAYATYATR\n", '-e', 'BamHI', '-a', 'possible', '-c')
        assert "no cut possible under any resolution" in out
        assert "pass --circular" not in out

    def test_gaps_reported_and_stripped(self):
        out = run_cli(">aln\nGG-AAT--ATTGG\n", '-e', 'SspI', '-a', 'definite')
        assert "10 bp, linear (3 gap character(s) removed)" in out
        assert "cut sites at: [5]" in out

    def test_state_does_not_leak_between_records(self):
        out = run_cli(">amb\nGGAAYATYATR\n>clean\nAAGAATTCTTGGATCCAA\n",
                      '-e', 'SspI,EcoRI', '-a', 'possible')
        _, clean = out.split("=== Digest of clean ===")
        assert "?" not in clean
        assert "maximal-cut scenario" not in clean
