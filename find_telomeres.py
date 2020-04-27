#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Reads a fasta file and uses regex to find telomere repeats.
'''

import argparse
import re

parser = argparse.ArgumentParser(description="Reads a fasta file and uses regex to find telomere repeats.")
parser.add_argument("reference", help="Reference genome fasta file. Required.", type = str)
args = parser.parse_args()

def readFasta():
    fa = {}
    with open(args.reference, "r") as fasta:
        sequence = None
        for line in fasta:
            line = line.strip()

            if line.startswith(">") and sequence == None:
                header = line[1:]
                sequence = []

            elif line.startswith(">") and sequence != None:
                # If new fasta entry, add old one to dict, pick new header and reset sequence
                fa[header] = "".join(sequence)
                header = line[1:]
                sequence = []

            else:
                sequence.append(line)

        # Last passthrough won't have any new entries, just add the remaining sequence
        fa[header] = "".join(sequence)

    return fa


def findTelomere(seq):
    '''
    Takes nucleotide sequence and checks if the sequence contains telomere repeats.
    '''

    # Look within first and last 200 bp for repeats
    start = seq[:100]
    end = seq[-100:]

    forward, reverse = False, False

    # Look for TAACCCC... and TTAGGGG...
    if re.search("TAA[C]+TAA[C]+", start):
        forward = True
    if re.search("TTA[G]+TTA[G]+", end):
        reverse = True

    return forward, reverse

def main():
    fasta = readFasta()
    n_f, n_r = 0,0

    for header, seq in fasta.items():
        forward, reverse = findTelomere(seq)

        if forward == True:
            print("{}\tforward\t{}".format(header, seq[:100]))
            n_f += 1
        if reverse == True:
            print("{}\treverse\t{}".format(header, seq[-100:]))
            n_r += 1

    print("\nTelomeres found: {} ({} forward, {} reverse)".format(str(n_f+n_r),n_f,n_r))

if __name__ == "__main__":
    main()
