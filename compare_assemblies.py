#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""compare_assemblies.py

Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Description: Calculate the number of overlapping bases between a reference and
    a query genome assembly.

Copyright (c) 2022, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
"""

__version__ = "0.1"

import mappy as mp
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-q", "--query_fasta", \
                    help="Query fasta file.", \
                    type = str)
parser.add_argument("-r", "--reference_fasta", \
                    help="Reference fasta file.", \
                    type = str)
parser.add_argument("-a","--alignment_preset", \
                    help="Minimap2 alignment present. Valid values: 'asm5', 'asm10', 'asm20'", \
                    default = "asm5", type = str)
parser.add_argument("-v","--version", help="Print version and quit.", action = "version", version = "compare_assemblies v.{}".format(__version__))
args = parser.parse_args()

def whole_genome_aligner(query,ref,alignment_preset):
    ''' Align two genomes and return the number of aligned bases
    '''
    idx = mp.Aligner(fn_idx_in=ref, preset=alignment_preset)
    n_bp = 0
    #import ipdb; ipdb.set_trace()

    for entry in query:
        seq = entry[1]
        for aln in idx.map(seq):
            if aln.is_primary:
                n_bp += aln.blen
    return n_bp

def printout(bp_q, bp_r, bp_ovl):
    ''' Print results to stdout.
    '''
    prcnt = round(bp_ovl/bp_r, 4)*100
    print("Reference: {} bp\nQuery: {} bp\nOverlapping: {} bp ({}% of reference)".format(bp_r, bp_q, bp_ovl, prcnt))


def main():
    query_fa = mp.fastx_read(args.query_fasta)
    ref_fa = mp.fastx_read(args.reference_fasta)
    n_bp_overlap = whole_genome_aligner(mp.fastx_read(args.query_fasta), args.reference_fasta, args.alignment_preset)
    n_bp_q = len( "".join([ entry[1] for entry in query_fa]) )
    n_bp_r = len( "".join([ entry[1] for entry in ref_fa]) )
    printout(n_bp_q, n_bp_r, n_bp_overlap)

if __name__ == "__main__":
    main()
