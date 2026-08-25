"""Static reference data for the restriction enzyme digest simulator.

Lookup tables only - no logic, and nothing imported from ``enzyme_digest``, so
this module can be read, diffed and regenerated on its own. ``enzyme_digest``
turns :data:`ENZYME_TABLE` into the typed :class:`~enzyme_digest.Enzyme`
records it works with.
"""

from typing import Dict, List, Tuple

# IUPAC nucleotide codes -> the set of bases each one stands for.
IUPAC_BASES: Dict[str, str] = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
IUPAC_ALPHABET = ''.join(IUPAC_BASES)

# Complementing a code complements the bases it stands for: R (A or G) pairs
# with Y (T or C), B (not A) with V (not T), and so on.
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')

# Recognition specificities taken from REBASE (rebase.neb.com) version 608, the
# database behind NEB's alphabetized list of recognition specificities.
#
# Values are (site, cut_top, cut_bottom): 0-based offsets from the first base of
# the site, both measured along the top strand, with the cut falling immediately
# before the base at that offset. See the Enzyme docstring in enzyme_digest.py.
# The comment on each line is the same enzyme in NEB notation.
ENZYME_TABLE: Dict[str, Tuple[str, int, int]] = {
    "EcoRI": ("GAATTC", 1, 5),              # G^AATTC
    "BamHI": ("GGATCC", 1, 5),              # G^GATCC
    "HindIII": ("AAGCTT", 1, 5),            # A^AGCTT
    "EcoRV": ("GATATC", 3, 3),              # GAT^ATC
    "PstI": ("CTGCAG", 5, 1),               # CTGCA^G
    "SmaI": ("CCCGGG", 3, 3),               # CCC^GGG
    "XbaI": ("TCTAGA", 1, 5),               # T^CTAGA
    "NotI": ("GCGGCCGC", 2, 6),             # GC^GGCCGC
    "SacI": ("GAGCTC", 5, 1),               # GAGCT^C
    "KpnI": ("GGTACC", 5, 1),               # GGTAC^C
    "SpeI": ("ACTAGT", 1, 5),               # A^CTAGT
    "BglII": ("AGATCT", 1, 5),              # A^GATCT
    "NcoI": ("CCATGG", 1, 5),               # C^CATGG
    "NdeI": ("CATATG", 2, 4),               # CA^TATG
    "XhoI": ("CTCGAG", 1, 5),               # C^TCGAG
    "SalI": ("GTCGAC", 1, 5),               # G^TCGAC
    "ClaI": ("ATCGAT", 2, 4),               # AT^CGAT
    "HpaI": ("GTTAAC", 3, 3),               # GTT^AAC
    "MluI": ("ACGCGT", 1, 5),               # A^CGCGT
    "NheI": ("GCTAGC", 1, 5),               # G^CTAGC
    "AatII": ("GACGTC", 5, 1),              # GACGT^C
    "AccI": ("GTMKAC", 2, 4),               # GT^MKAC
    "AflII": ("CTTAAG", 1, 5),              # C^TTAAG
    "AgeI": ("ACCGGT", 1, 5),               # A^CCGGT
    "ApaI": ("GGGCCC", 5, 1),               # GGGCC^C
    "ApaLI": ("GTGCAC", 1, 5),              # G^TGCAC
    "AscI": ("GGCGCGCC", 2, 6),             # GG^CGCGCC
    "AvaI": ("CYCGRG", 1, 5),               # C^YCGRG
    "AvrII": ("CCTAGG", 1, 5),              # C^CTAGG
    "BclI": ("TGATCA", 1, 5),               # T^GATCA
    "BsrGI": ("TGTACA", 1, 5),              # T^GTACA
    "BstEII": ("GGTNACC", 1, 6),            # G^GTNACC
    "BstXI": ("CCANNNNNNTGG", 8, 4),        # CCANNNNN^NTGG
    "DraI": ("TTTAAA", 3, 3),               # TTT^AAA
    "DraIII": ("CACNNNGTG", 6, 3),          # CACNNN^GTG
    "FseI": ("GGCCGGCC", 6, 2),             # GGCCGG^CC
    "HaeIII": ("GGCC", 2, 2),               # GG^CC
    "HincII": ("GTYRAC", 3, 3),             # GTY^RAC
    "HinfI": ("GANTC", 1, 4),               # G^ANTC
    "HpaII": ("CCGG", 1, 3),                # C^CGG
    "MspI": ("CCGG", 1, 3),                 # C^CGG
    "MseI": ("TTAA", 1, 3),                 # T^TAA
    "MfeI": ("CAATTG", 1, 5),               # C^AATTG
    "NarI": ("GGCGCC", 2, 4),               # GG^CGCC
    "NruI": ("TCGCGA", 3, 3),               # TCG^CGA
    "PacI": ("TTAATTAA", 5, 3),             # TTAAT^TAA
    "PmeI": ("GTTTAAAC", 4, 4),             # GTTT^AAAC
    "PvuI": ("CGATCG", 4, 2),               # CGAT^CG
    "PvuII": ("CAGCTG", 3, 3),              # CAG^CTG
    "RsaI": ("GTAC", 2, 2),                 # GT^AC
    "SacII": ("CCGCGG", 4, 2),              # CCGC^GG
    "ScaI": ("AGTACT", 3, 3),               # AGT^ACT
    "SfiI": ("GGCCNNNNNGGCC", 8, 5),        # GGCCNNNN^NGGCC
    "SgrAI": ("CRCCGGYG", 2, 6),            # CR^CCGGYG
    "SnaBI": ("TACGTA", 3, 3),              # TAC^GTA
    "SphI": ("GCATGC", 5, 1),               # GCATG^C
    "SspI": ("AATATT", 3, 3),               # AAT^ATT
    "StuI": ("AGGCCT", 3, 3),               # AGG^CCT
    "StyI": ("CCWWGG", 1, 5),               # C^CWWGG
    "SwaI": ("ATTTAAAT", 4, 4),             # ATTT^AAAT
    "TaqI": ("TCGA", 1, 3),                 # T^CGA
    "XcmI": ("CCANNNNNNNNNTGG", 8, 7),      # CCANNNNN^NNNNTGG
    "XmaI": ("CCCGGG", 1, 5),               # C^CCGGG
    "ZraI": ("GACGTC", 3, 3),               # GAC^GTC
    "AlwNI": ("CAGNNNCTG", 6, 3),           # CAGNNN^CTG
    "AhdI": ("GACNNNNNGTC", 6, 5),          # GACNNN^NNGTC
    "XmnI": ("GAANNNNTTC", 5, 5),           # GAANN^NNTTC
    "ApoI": ("RAATTY", 1, 5),               # R^AATTY
    "NspI": ("RCATGY", 5, 1),               # RCATG^Y
    "HaeII": ("RGCGCY", 5, 1),              # RGCGC^Y
    "MboI": ("GATC", 0, 4),                 # ^GATC
    "Sau3AI": ("GATC", 0, 4),               # ^GATC
    "BsaI": ("GGTCTC", 7, 11),              # GGTCTC(1/5)
    "BsmBI": ("CGTCTC", 7, 11),             # CGTCTC(1/5)
    "BbsI": ("GAAGAC", 8, 12),              # GAAGAC(2/6)
    "SapI": ("GCTCTTC", 8, 11),             # GCTCTTC(1/4)
    "BspQI": ("GCTCTTC", 8, 11),            # GCTCTTC(1/4)
    "BsmAI": ("GTCTC", 6, 10),              # GTCTC(1/5)
    "AlwI": ("GGATC", 9, 10),               # GGATC(4/5)
    "PleI": ("GAGTC", 9, 10),               # GAGTC(4/5)
    "MlyI": ("GAGTC", 10, 10),              # GAGTC(5/5)
    "HgaI": ("GACGC", 10, 15),              # GACGC(5/10)
    "BbvI": ("GCAGC", 13, 17),              # GCAGC(8/12)
    "FokI": ("GGATG", 14, 18),              # GGATG(9/13)
    "BtgZI": ("GCGATG", 16, 20),            # GCGATG(10/14)
    "AarI": ("CACCTGC", 11, 15),            # CACCTGC(4/8)
    "BciVI": ("GTATCC", 12, 11),            # GTATCC(6/5)
    "BsgI": ("GTGCAG", 22, 20),             # GTGCAG(16/14)
}

# Band sizes of a standard 100 bp DNA ladder, for the simulated gel.
DNA_LADDER_100BP: List[int] = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1200, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000
]
