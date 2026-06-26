import pytest
import tempfile
import os
from enzyme_digest import (
    parse_fasta,
    parse_enzyme_spec,
    find_cut_positions,
    digest_linear,
    digest_circular,
    ENZYME_DB
)


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
        name, recog, off = parse_enzyme_spec("EcoRI")
        assert name == "ECORI"
        assert recog == "GAATTC"
        assert off == 1

    def test_case_insensitive_builtin(self):
        name, recog, off = parse_enzyme_spec("ecori")
        assert name == "ECORI"

    def test_custom_valid(self):
        name, recog, off = parse_enzyme_spec("MyI:ACG:2")
        assert name == "MYI"
        assert recog == "ACG"
        assert off == 2

    def test_custom_missing_offset_raises(self):
        with pytest.raises(ValueError, match="offset required"):
            parse_enzyme_spec("Test:GAATTC")

    def test_custom_offset_out_of_range_raises(self):
        with pytest.raises(ValueError, match="out of range"):
            parse_enzyme_spec("Bad:ACG:5")

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError, match="Unknown enzyme"):
            parse_enzyme_spec("NonExistent")


class TestFindCutPositions:
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
