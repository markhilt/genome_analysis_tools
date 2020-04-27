#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
convert_mol_cov.py
Version 0.2
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

This script reads a bed file as output by tigmint to sum of the molecule
coverage per base pair.

This version is a bug fix from the paper version.

Copyright (c) 2019, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Reads a bed file as output by tigmint to infer molecule sizes.")
parser.add_argument("bed", help="Bed file. Required.", type = str)
parser.add_argument("faidx", help="Fasta index file. Required.", type = str)
parser.add_argument("-s", "--stepsize", help="Interval between positions to check coverage", type = int)
#parser.add_argument("-m", "--mode", help="Mode: sort or ", type = str)

args = parser.parse_args()

def reportProgress(current,total):
    return "Completed: {0}% ({1} out of {2})".format( str(round( (current / total) * 100, 2)), current, total)


def passBed():
    '''
    Goes through a bed file to collect general information about contigs and sizes
    '''
    asm = {} # {tig: <len of tig>}
    with open(args.bed, "r") as bed:
        for line in bed:
            line = line.strip()

            tig = line.split("\t")[0]
            coord = int(line.split("\t")[2])

            if tig not in asm.keys():
                asm[tig] = coord

            elif coord > asm[tig]:
                asm[tig] = coord

    return asm


def readFaidx(faidx):
    ''' Alternative for genome stats.
    '''
    asm = {} # {tig: <len of tig>}

    with open(faidx, "r") as idx:
        for line in idx:
            line = line.strip()
            tig, len = line.split("\t")[0], int(line.split("\t")[1])
            asm[tig] = len
    return asm


def bed2cov(genome):
    '''
    Goes through a bed file to sum up the coverage at every base position
    '''
    cov = {} # This will be a dict of numpy arrays
    for tig in genome.keys():
        cov[tig] = np.zeros(genome[tig])
    with open(args.bed, "r") as bed:
        # Fill the array with the physical coverage from the bed file
        for line in bed:
            tig = line.split("\t")[0]
            molspan = (int(line.split("\t")[1]), int(line.split("\t")[2]))
            cov[tig][molspan[0]-1:molspan[1]] += 1 # Increment coverage at these sites by one
    return cov

def main():
    print("Collecting genome stats...")
    genome = readFaidx(args.faidx)
    print("Done\n\nReading molecule coverage...")
    coverage = bed2cov(genome)
    print("Done\n\nWriting output...")

    with open(args.bed.split(".bed")[0] + ".mol.cov", "w") as out:
        for k,v in coverage.items():
            c = 1
            for val in v:
                out.write("{}\t{}\t{}\n".format(k,str(c),val))
                c += 1

if __name__ == "__main__":
    main()
