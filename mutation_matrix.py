#!/usr/bin/python
#


# The perl script is used only once. no need to rewrite it in python.

import argparse
from Bio import SeqIO


usage="""

Generate a matrix of mutation counts between all possible codons. Each column/row of the matrix corresponds to one codon.
The matrix entries indicate how many mutations are at least necessary to get from one codon to the other.
The matrix is generated once and then used as reference by todb_mutations_from_align.

The underlying model for mutations assumes the lowest number of mutations and does not differentiate between nucleotides.

AUTHOR /n Francisco Arcila rewriten from Katharina Imkeller

"""




parser = argparse.ArgumentParser(description=usage)
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action="store_true")
#group.add_argument("-q", "--quiet", action="store_true")
#parser.add_argument("x", type=int, help="the base")
#parser.add_argument("y", type=int, help="the exponent")
args = parser.parse_args()
#answer = args.x**args.y


def hamdist (str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs


# dna1 = 0
# dna2 = 0



# generate matrix with all codons
characters = ["A", "C", "G", "T", "-"]
codons = []
poss_parents = []
poss_daughters = []


count = 0
# create a list with all possible nucleotides
for nt1 in characters:
    for nt2 in characters:
        for nt3 in characters:
            codons.append(nt1+nt2+nt3)

print (codons)
for obs_codon in codons:
    for ger_codon in codons:
        insertion, deletion, mutation, n_repl, n_silent = 0,0,0,0,0

        weight = []
