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

# --- BEGIN GENERATED TABLES (tools/build_enzyme_table.py) ---
# Generated from restriction_enzymes.csv - do not edit by hand.
# Regenerate with: python tools/build_enzyme_table.py
#
# (site, cut_top, cut_bottom): 0-based offsets from the first base of the
# site, both measured along the top strand, the cut falling immediately
# before the base at that offset. Offsets may sit outside the site - type
# IIS enzymes cut downstream of it. The comment is the NEB specificity.
ENZYME_TABLE: Dict[str, Tuple[str, int, int]] = {
    "AatII":     ('GACGTC', 5, 1),                             # GACGT^C
    "AbsI":      ('CCTCGAGG', 2, 6),                           # CC^TCGAGG
    "Acc65I":    ('GGTACC', 1, 5),                             # G^GTACC
    "AccI":      ('GTMKAC', 2, 4),                             # GT^MKAC
    "AciI":      ('CCGC', 1, 3),                               # C^CGC
    "AclI":      ('AACGTT', 2, 4),                             # AA^CGTT
    "AcuI":      ('CTGAAG', 22, 20),                           # CTGAAG(16/14)
    "AfeI":      ('AGCGCT', 3, 3),                             # AGC^GCT
    "AflII":     ('CTTAAG', 1, 5),                             # C^TTAAG
    "AflIII":    ('ACRYGT', 1, 5),                             # A^CRYGT
    "AgeI":      ('ACCGGT', 1, 5),                             # A^CCGGT
    "AgsI":      ('TTSAA', 3, 2),                              # TTS^AA
    "AhdI":      ('GACNNNNNGTC', 6, 5),                        # GACNNN^NNGTC
    "AleI":      ('CACNNNNGTG', 5, 5),                         # CACNN^NNGTG
    "AluI":      ('AGCT', 2, 2),                               # AG^CT
    "AlwI":      ('GGATC', 9, 10),                             # GGATC(4/5)
    "AlwNI":     ('CAGNNNCTG', 6, 3),                          # CAGNNN^CTG
    "AoxI":      ('GGCC', 0, 4),                               # ^GGCC
    "ApaI":      ('GGGCCC', 5, 1),                             # GGGCC^C
    "ApaLI":     ('GTGCAC', 1, 5),                             # G^TGCAC
    "ApeKI":     ('GCWGC', 1, 4),                              # G^CWGC
    "ApoI":      ('RAATTY', 1, 5),                             # R^AATTY
    "AscI":      ('GGCGCGCC', 2, 6),                           # GG^CGCGCC
    "AseI":      ('ATTAAT', 2, 4),                             # AT^TAAT
    "AsiSI":     ('GCGATCGC', 5, 3),                           # GCGAT^CGC
    "AvaI":      ('CYCGRG', 1, 5),                             # C^YCGRG
    "AvaII":     ('GGWCC', 1, 4),                              # G^GWCC
    "AvrII":     ('CCTAGG', 1, 5),                             # C^CTAGG
    "BaeGI":     ('GKGCMC', 5, 1),                             # GKGCM^C
    "BamHI":     ('GGATCC', 1, 5),                             # G^GATCC
    "BanI":      ('GGYRCC', 1, 5),                             # G^GYRCC
    "BanII":     ('GRGCYC', 5, 1),                             # GRGCY^C
    "BbsI":      ('GAAGAC', 8, 12),                            # GAAGAC(2/6)
    "BbvCI":     ('CCTCAGC', 2, 5),                            # CC^TCAGC
    "BbvI":      ('GCAGC', 13, 17),                            # GCAGC(8/12)
    "BccI":      ('CCATC', 9, 10),                             # CCATC(4/5)
    "BceAI":     ('ACGGC', 17, 19),                            # ACGGC(12/14)
    "BciVI":     ('GTATCC', 12, 11),                           # GTATCC(6/5)
    "BclI":      ('TGATCA', 1, 5),                             # T^GATCA
    "BcoDI":     ('GTCTC', 6, 10),                             # GTCTC(1/5)
    "BfaI":      ('CTAG', 1, 3),                               # C^TAG
    "BfuAI":     ('ACCTGC', 10, 14),                           # ACCTGC(4/8)
    "BglI":      ('GCCNNNNNGGC', 7, 4),                        # GCCNNNN^NGGC
    "BglII":     ('AGATCT', 1, 5),                             # A^GATCT
    "BlpI":      ('GCTNAGC', 2, 5),                            # GC^TNAGC
    "BlsI":      ('GCNGC', 3, 2),                              # GCN^GC
    "BmgBI":     ('CACGTC', 3, 3),                             # CAC^GTC
    "BmrI":      ('ACTGGG', 11, 10),                           # ACTGGG(5/4)
    "BmtI":      ('GCTAGC', 5, 1),                             # GCTAG^C
    "BpmI":      ('CTGGAG', 22, 20),                           # CTGGAG(16/14)
    "Bpu10I":    ('CCTNAGC', 2, 5),                            # CC^TNAGC
    "BpuEI":     ('CTTGAG', 22, 20),                           # CTTGAG(16/14)
    "BsaAI":     ('YACGTR', 3, 3),                             # YAC^GTR
    "BsaBI":     ('GATNNNNATC', 5, 5),                         # GATNN^NNATC
    "BsaHI":     ('GRCGYC', 2, 4),                             # GR^CGYC
    "BsaI":      ('GGTCTC', 7, 11),                            # GGTCTC(1/5)
    "BsaJI":     ('CCNNGG', 1, 5),                             # C^CNNGG
    "BsaWI":     ('WCCGGW', 1, 5),                             # W^CCGGW
    "BseMII":    ('CTCAG', 15, 13),                            # CTCAG(10/8)
    "BseRI":     ('GAGGAG', 16, 14),                           # GAGGAG(10/8)
    "BseYI":     ('CCCAGC', 1, 5),                             # C^CCAGC
    "BsgI":      ('GTGCAG', 22, 20),                           # GTGCAG(16/14)
    "BsiEI":     ('CGRYCG', 4, 2),                             # CGRY^CG
    "BsiHKAI":   ('GWGCWC', 5, 1),                             # GWGCW^C
    "BsiWI":     ('CGTACG', 1, 5),                             # C^GTACG
    "BslI":      ('CCNNNNNNNGG', 7, 4),                        # CCNNNNN^NNGG
    "BsmAI":     ('GTCTC', 6, 10),                             # GTCTC(1/5)
    "BsmBI":     ('CGTCTC', 7, 11),                            # CGTCTC(1/5)
    "BsmFI":     ('GGGAC', 15, 19),                            # GGGAC(10/14)
    "BsmI":      ('GAATGC', 7, 5),                             # GAATGC(1/-1)
    "BsoBI":     ('CYCGRG', 1, 5),                             # C^YCGRG
    "Bsp1286I":  ('GDGCHC', 5, 1),                             # GDGCH^C
    "BspCNI":    ('CTCAG', 14, 12),                            # CTCAG(9/7)
    "BspDI":     ('ATCGAT', 2, 4),                             # AT^CGAT
    "BspEI":     ('TCCGGA', 1, 5),                             # T^CCGGA
    "BspHI":     ('TCATGA', 1, 5),                             # T^CATGA
    "BspMI":     ('ACCTGC', 10, 14),                           # ACCTGC(4/8)
    "BspQI":     ('GCTCTTC', 8, 11),                           # GCTCTTC(1/4)
    "BsrBI":     ('CCGCTC', 3, 3),                             # CCG^CTC
    "BsrDI":     ('GCAATG', 8, 6),                             # GCAATG(2/0)
    "BsrFI":     ('RCCGGY', 1, 5),                             # R^CCGGY
    "BsrGI":     ('TGTACA', 1, 5),                             # T^GTACA
    "BsrI":      ('ACTGG', 6, 4),                              # ACTGG(1/-1)
    "BssHII":    ('GCGCGC', 1, 5),                             # G^CGCGC
    "BssSI":     ('CACGAG', 1, 5),                             # C^ACGAG
    "BstAPI":    ('GCANNNNNTGC', 7, 4),                        # GCANNNN^NTGC
    "BstBI":     ('TTCGAA', 2, 4),                             # TT^CGAA
    "BstEII":    ('GGTNACC', 1, 6),                            # G^GTNACC
    "BstKTI":    ('GATC', 3, 1),                               # GAT^C
    "BstNI":     ('CCWGG', 2, 3),                              # CC^WGG
    "BstUI":     ('CGCG', 2, 2),                               # CG^CG
    "BstXI":     ('CCANNNNNNTGG', 8, 4),                       # CCANNNNN^NTGG
    "BstYI":     ('RGATCY', 1, 5),                             # R^GATCY
    "BstZ17I":   ('GTATAC', 3, 3),                             # GTA^TAC
    "Bsu36I":    ('CCTNAGG', 2, 5),                            # CC^TNAGG
    "BtgI":      ('CCRYGG', 1, 5),                             # C^CRYGG
    "BtgZI":     ('GCGATG', 16, 20),                           # GCGATG(10/14)
    "BtsCI":     ('GGATG', 7, 5),                              # GGATG(2/0)
    "BtsI":      ('GCAGTG', 8, 6),                             # GCAGTG(2/0)
    "BtsIMutI":  ('CAGTG', 7, 5),                              # CAGTG(2/0)
    "Cac8I":     ('GCNNGC', 3, 3),                             # GCN^NGC
    "ClaI":      ('ATCGAT', 2, 4),                             # AT^CGAT
    "CviKI-1":   ('RGCY', 2, 2),                               # RG^CY
    "CviQI":     ('GTAC', 1, 3),                               # G^TAC
    "DdeI":      ('CTNAG', 1, 4),                              # C^TNAG
    "DpnI":      ('GATC', 2, 2),                               # GA^TC
    "DpnII":     ('GATC', 0, 4),                               # ^GATC
    "DraI":      ('TTTAAA', 3, 3),                             # TTT^AAA
    "DraIII":    ('CACNNNGTG', 6, 3),                          # CACNNN^GTG
    "DrdI":      ('GACNNNNNNGTC', 7, 5),                       # GACNNNN^NNGTC
    "EaeI":      ('YGGCCR', 1, 5),                             # Y^GGCCR
    "EagI":      ('CGGCCG', 1, 5),                             # C^GGCCG
    "EarI":      ('CTCTTC', 7, 10),                            # CTCTTC(1/4)
    "EciI":      ('GGCGGA', 17, 15),                           # GGCGGA(11/9)
    "Eco53kI":   ('GAGCTC', 3, 3),                             # GAG^CTC
    "EcoNI":     ('CCTNNNNNAGG', 5, 6),                        # CCTNN^NNNAGG
    "EcoO109I":  ('RGGNCCY', 2, 5),                            # RG^GNCCY
    "EcoP15I":   ('CAGCAG', 31, 33),                           # CAGCAG(25/27)
    "EcoRI":     ('GAATTC', 1, 5),                             # G^AATTC
    "EcoRV":     ('GATATC', 3, 3),                             # GAT^ATC
    "Esp3I":     ('CGTCTC', 7, 11),                            # CGTCTC(1/5)
    "FaiI":      ('YATR', 2, 2),                               # YA^TR
    "FatI":      ('CATG', 0, 4),                               # ^CATG
    "FauI":      ('CCCGC', 9, 11),                             # CCCGC(4/6)
    "Fnu4HI":    ('GCNGC', 2, 3),                              # GC^NGC
    "FokI":      ('GGATG', 14, 18),                            # GGATG(9/13)
    "FseI":      ('GGCCGGCC', 6, 2),                           # GGCCGG^CC
    "FspAI":     ('RTGCGCAY', 4, 4),                           # RTGC^GCAY
    "FspI":      ('TGCGCA', 3, 3),                             # TGC^GCA
    "GlaI":      ('GCGC', 2, 2),                               # GC^GC
    "GsaI":      ('CCCAGC', 5, 1),                             # CCCAG^C
    "HaeII":     ('RGCGCY', 5, 1),                             # RGCGC^Y
    "HaeIII":    ('GGCC', 2, 2),                               # GG^CC
    "HgaI":      ('GACGC', 10, 15),                            # GACGC(5/10)
    "HhaI":      ('GCGC', 3, 1),                               # GCG^C
    "HincII":    ('GTYRAC', 3, 3),                             # GTY^RAC
    "HindIII":   ('AAGCTT', 1, 5),                             # A^AGCTT
    "HinfI":     ('GANTC', 1, 4),                              # G^ANTC
    "HinP1I":    ('GCGC', 1, 3),                               # G^CGC
    "HpaI":      ('GTTAAC', 3, 3),                             # GTT^AAC
    "HpaII":     ('CCGG', 1, 3),                               # C^CGG
    "HphI":      ('GGTGA', 13, 12),                            # GGTGA(8/7)
    "Hpy166II":  ('GTNNAC', 3, 3),                             # GTN^NAC
    "Hpy188I":   ('TCNGA', 3, 2),                              # TCN^GA
    "Hpy188III": ('TCNNGA', 2, 4),                             # TC^NNGA
    "Hpy99I":    ('CGWCG', 5, 0),                              # CGWCG^
    "HpyAV":     ('CCTTC', 11, 10),                            # CCTTC(6/5)
    "HpyCH4III": ('ACNGT', 3, 2),                              # ACN^GT
    "HpyCH4IV":  ('ACGT', 1, 3),                               # A^CGT
    "HpyCH4V":   ('TGCA', 2, 2),                               # TG^CA
    "I-CeuI":    ('TAACTATAACGGTCCTAAGGTAGCGAA', 18, 14),      # TAACTATAACGGTCCTAAGGTAGCGAA(-9/-13)
    "I-PpoI":    ('TAACTATGACTCTCTTAAGGTAGCCAAAT', 18, 14),    # TAACTATGACTCTCTTAAGGTAGCCAAAT(-11/-15)
    "I-SceI":    ('TAGGGATAACAGGGTAAT', 9, 5),                 # TAGGGATAACAGGGTAAT(-9/-13)
    "KasI":      ('GGCGCC', 1, 5),                             # G^GCGCC
    "KpnI":      ('GGTACC', 5, 1),                             # GGTAC^C
    "LmnI":      ('GCTCC', 6, 4),                              # GCTCC(1/-1)
    "MauBI":     ('CGCGCGCG', 2, 6),                           # CG^CGCGCG
    "MboI":      ('GATC', 0, 4),                               # ^GATC
    "MboII":     ('GAAGA', 13, 12),                            # GAAGA(8/7)
    "MfeI":      ('CAATTG', 1, 5),                             # C^AATTG
    "MluCI":     ('AATT', 0, 4),                               # ^AATT
    "MluI":      ('ACGCGT', 1, 5),                             # A^CGCGT
    "MlyI":      ('GAGTC', 10, 10),                            # GAGTC(5/5)
    "MmeI":      ('TCCRAC', 26, 24),                           # TCCRAC(20/18)
    "MnlI":      ('CCTC', 11, 10),                             # CCTC(7/6)
    "MreI":      ('CGCCGGCG', 2, 6),                           # CG^CCGGCG
    "MscI":      ('TGGCCA', 3, 3),                             # TGG^CCA
    "MseI":      ('TTAA', 1, 3),                               # T^TAA
    "MslI":      ('CAYNNNNRTG', 5, 5),                         # CAYNN^NNRTG
    "MspA1I":    ('CMGCKG', 3, 3),                             # CMG^CKG
    "MspI":      ('CCGG', 1, 3),                               # C^CGG
    "MspJI":     ('CNNR', 13, 17),                             # CNNR(9/13)
    "MteI":      ('GCGCNGCGC', 4, 5),                          # GCGC^NGCGC
    "MwoI":      ('GCNNNNNNNGC', 7, 4),                        # GCNNNNN^NNGC
    "NaeI":      ('GCCGGC', 3, 3),                             # GCC^GGC
    "NarI":      ('GGCGCC', 2, 4),                             # GG^CGCC
    "NciI":      ('CCSGG', 2, 3),                              # CC^SGG
    "NcoI":      ('CCATGG', 1, 5),                             # C^CATGG
    "NdeI":      ('CATATG', 2, 4),                             # CA^TATG
    "NgoMIV":    ('GCCGGC', 1, 5),                             # G^CCGGC
    "NheI":      ('GCTAGC', 1, 5),                             # G^CTAGC
    "NlaIII":    ('CATG', 4, 0),                               # CATG^
    "NlaIV":     ('GGNNCC', 3, 3),                             # GGN^NCC
    "NmeAIII":   ('GCCGAG', 27, 25),                           # GCCGAG(21/19)
    "NotI":      ('GCGGCCGC', 2, 6),                           # GC^GGCCGC
    "NruI":      ('TCGCGA', 3, 3),                             # TCG^CGA
    "NsiI":      ('ATGCAT', 5, 1),                             # ATGCA^T
    "NspI":      ('RCATGY', 5, 1),                             # RCATG^Y
    "PacI":      ('TTAATTAA', 5, 3),                           # TTAAT^TAA
    "PaeR7I":    ('CTCGAG', 1, 5),                             # C^TCGAG
    "PaqCI":     ('CACCTGC', 11, 15),                          # CACCTGC(4/8)
    "PasI":      ('CCCWGGG', 2, 5),                            # CC^CWGGG
    "PciI":      ('ACATGT', 1, 5),                             # A^CATGT
    "PcsI":      ('WCGNNNNNNNCGW', 7, 6),                      # WCGNNNN^NNNCGW
    "PflFI":     ('GACNNNGTC', 4, 5),                          # GACN^NNGTC
    "PflMI":     ('CCANNNNNTGG', 7, 4),                        # CCANNNN^NTGG
    "PfoI":      ('TCCNGGA', 1, 6),                            # T^CCNGGA
    "PI-PspI":   ('TGGCAAACAGCTATTATGGGTATTATGGGT', 17, 13),   # TGGCAAACAGCTATTAT^GGGTATTATGGGT
    "PI-SceI":   ('ATCTATGTCGGGTGCGGAGAAAGAGGTAATGAAATGG', 15, 11), # ATCTATGTCGGGTGCGGAGAAAGAGGTAATGAAATGG(-22/-26)
    "PleI":      ('GAGTC', 9, 10),                             # GAGTC(4/5)
    "PluTI":     ('GGCGCC', 5, 1),                             # GGCGC^C
    "PmeI":      ('GTTTAAAC', 4, 4),                           # GTTT^AAAC
    "PmlI":      ('CACGTG', 3, 3),                             # CAC^GTG
    "PpuMI":     ('RGGWCCY', 2, 5),                            # RG^GWCCY
    "PshAI":     ('GACNNNNGTC', 5, 5),                         # GACNN^NNGTC
    "PsiI":      ('TTATAA', 3, 3),                             # TTA^TAA
    "PspGI":     ('CCWGG', 0, 5),                              # ^CCWGG
    "PspOMI":    ('GGGCCC', 1, 5),                             # G^GGCCC
    "PspXI":     ('VCTCGAGB', 2, 6),                           # VC^TCGAGB
    "PstI":      ('CTGCAG', 5, 1),                             # CTGCA^G
    "PvuI":      ('CGATCG', 4, 2),                             # CGAT^CG
    "PvuII":     ('CAGCTG', 3, 3),                             # CAG^CTG
    "RsaI":      ('GTAC', 2, 2),                               # GT^AC
    "RsrII":     ('CGGWCCG', 2, 5),                            # CG^GWCCG
    "SacI":      ('GAGCTC', 5, 1),                             # GAGCT^C
    "SacII":     ('CCGCGG', 4, 2),                             # CCGC^GG
    "SalI":      ('GTCGAC', 1, 5),                             # G^TCGAC
    "SapI":      ('GCTCTTC', 8, 11),                           # GCTCTTC(1/4)
    "Sau3AI":    ('GATC', 0, 4),                               # ^GATC
    "Sau96I":    ('GGNCC', 1, 4),                              # G^GNCC
    "SbfI":      ('CCTGCAGG', 6, 2),                           # CCTGCA^GG
    "ScaI":      ('AGTACT', 3, 3),                             # AGT^ACT
    "ScrFI":     ('CCNGG', 2, 3),                              # CC^NGG
    "SetI":      ('ASST', 4, 0),                               # ASST^
    "SexAI":     ('ACCWGGT', 1, 6),                            # A^CCWGGT
    "SfaNI":     ('GCATC', 10, 14),                            # GCATC(5/9)
    "SfcI":      ('CTRYAG', 1, 5),                             # C^TRYAG
    "SfiI":      ('GGCCNNNNNGGCC', 8, 5),                      # GGCCNNNN^NGGCC
    "SfoI":      ('GGCGCC', 3, 3),                             # GGC^GCC
    "SgeI":      ('CNNG', 13, 17),                             # CNNG(9/13)
    "SgrAI":     ('CRCCGGYG', 2, 6),                           # CR^CCGGYG
    "SgrDI":     ('CGTCGACG', 2, 6),                           # CG^TCGACG
    "SmaI":      ('CCCGGG', 3, 3),                             # CCC^GGG
    "SmlI":      ('CTYRAG', 1, 5),                             # C^TYRAG
    "SnaBI":     ('TACGTA', 3, 3),                             # TAC^GTA
    "SpeI":      ('ACTAGT', 1, 5),                             # A^CTAGT
    "SphI":      ('GCATGC', 5, 1),                             # GCATG^C
    "SrfI":      ('GCCCGGGC', 4, 4),                           # GCCC^GGGC
    "SspI":      ('AATATT', 3, 3),                             # AAT^ATT
    "StuI":      ('AGGCCT', 3, 3),                             # AGG^CCT
    "StyD4I":    ('CCNGG', 0, 5),                              # ^CCNGG
    "StyI":      ('CCWWGG', 1, 5),                             # C^CWWGG
    "SwaI":      ('ATTTAAAT', 4, 4),                           # ATTT^AAAT
    "TaiI":      ('ACGT', 4, 0),                               # ACGT^
    "TaqI":      ('TCGA', 1, 3),                               # T^CGA
    "TaqII":     ('GACCGA', 17, 15),                           # GACCGA(11/9)
    "TatI":      ('WGTACW', 1, 5),                             # W^GTACW
    "TauI":      ('GCSGC', 4, 1),                              # GCSG^C
    "TfiI":      ('GAWTC', 1, 4),                              # G^AWTC
    "TseI":      ('GCWGC', 1, 4),                              # G^CWGC
    "Tsp45I":    ('GTSAC', 0, 5),                              # ^GTSAC
    "TspDTI":    ('ATGAA', 16, 14),                            # ATGAA(11/9)
    "TspGWI":    ('ACGGA', 16, 14),                            # ACGGA(11/9)
    "TspMI":     ('CCCGGG', 1, 5),                             # C^CCGGG
    "TspRI":     ('CASTG', 7, -2),                             # CASTG(2/-7)
    "Tth111I":   ('GACNNNGTC', 4, 5),                          # GACN^NNGTC
    "XbaI":      ('TCTAGA', 1, 5),                             # T^CTAGA
    "XcmI":      ('CCANNNNNNNNNTGG', 8, 7),                    # CCANNNNN^NNNNTGG
    "XhoI":      ('CTCGAG', 1, 5),                             # C^TCGAG
    "XmaI":      ('CCCGGG', 1, 5),                             # C^CCGGG
    "XmnI":      ('GAANNNNTTC', 5, 5),                         # GAANN^NNTTC
    "ZraI":      ('GACGTC', 3, 3),                             # GAC^GTC
}

# Nicking enzymes cut one strand only. They leave the duplex intact, so
# there is no fragment pattern to simulate; named here to give a better
# error than 'unknown enzyme'.
NICKING_ENZYMES: Dict[str, str] = {
    "Nb.BbvCI": "CCTCA_GC",
    "Nb.Bpu10I": "CCTNA_GC",
    "Nb.BsmI": "GAATG_C",
    "Nb.BsrDI": "GCAATG_",
    "Nb.BssSI": "CACGA_G",
    "Nb.BtsI": "GCAGTG_",
    "Nt.AlwI": "GGATCNNNN^",
    "Nt.BbvCI": "CC^TCAGC",
    "Nt.Bpu10I": "CC^TNAGC",
    "Nt.BsmAI": "GTCTCN^",
    "Nt.BspQI": "GCTCTTCN^",
    "Nt.BstNBI": "GAGTCNNNN^",
    "Nt.CviPII": "^CCD",
}

# Type IIB enzymes cut on both sides of their site, excising it. Two cuts
# per site do not fit the single-cut Enzyme record, so they are listed
# rather than loaded.
DUAL_CUT_ENZYMES: Dict[str, str] = {
    "AjuI": "_NNNNN^NNNNNNNGAANNNNNNNTTGGNNNNNN_NNNNN^",
    "AlfI": "_NN^NNNNNNNNNNGCANNNNNNTGCNNNNNNNNNN_NN^",
    "AloI": "_NNNNN^NNNNNNNGAACNNNNNNTCCNNNNNNN_NNNNN^",
    "ArsI": "_NNNNN^NNNNNNNNGACNNNNNNTTYGNNNNNN_NNNNN^",
    "BaeI": "_NNNNN^NNNNNNNNNNACNNNNGTAYCNNNNNNN_NNNNN^",
    "BarI": "_NNNNN^NNNNNNNGAAGNNNNNNTACNNNNNNN_NNNNN^",
    "BcgI": "_NN^NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNN_NN^",
    "BplI": "_NNNNN^NNNNNNNNGAGNNNNNCTCNNNNNNNN_NNNNN^",
    "BsaXI": "_NNN^NNNNNNNNNACNNNNNCTCCNNNNNNN_NNN^",
    "CspCI": "_NN^NNNNNNNNNNNCAANNNNNGTGGNNNNNNNNNN_NN^",
    "FalI": "_NNNNN^NNNNNNNNAAGNNNNNCTTNNNNNNNN_NNNNN^",
    "PsrI": "_NNNNN^NNNNNNNGAACNNNNNNTACNNNNNNN_NNNNN^",
}
# --- END GENERATED TABLES ---

# Band sizes of a standard 100 bp DNA ladder, for the simulated gel.
DNA_LADDER_100BP: List[int] = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1200, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000
]
