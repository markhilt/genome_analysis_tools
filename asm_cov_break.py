#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
convert_mol_cov.py
Version 0.3
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

This script reads a bed file as with coverage per base pair and finds regions of
low coverage, to split a given genome assembly.

Works in combination with output from convert_mol_cov.py, or any other coverage
bed files.

Copyright (c) 2020, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""


import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Reads a bed file to split an \
                                            assembly in low-covered regions.")
parser.add_argument("bed", help="Bed file. Required.", type = str)
parser.add_argument("fasta", help="Fasta file. Required.", type = str)
parser.add_argument("-c", "--coverage", help="Minimum coverage [10]", default = 10, type = int)
parser.add_argument("-o", "--out", help="Output prefix", default = "broken", type = str)
parser.add_argument("-l", "--length", help="Minimum contig length to consider for breaking [10000]", default = 10000, type = int)
parser.add_argument("-d", "--end_distance", help="Ignore regions up to this distance from contig ends [1000]", default = 1000, type = int)
args = parser.parse_args()

def reportProgress(current,total):
    return "Completed: {0}% ({1} out of {2})".format( str(round( (current / total) * 100, 2)), current, total)

def readCov(bed_cov, min_cov, min_len, distance):
    breaks = {}
    with open(bed_cov, "r") as bed:
        prev_tig = ""
        break_start, tig_start = 0, 0
        for line in bed:
            line = line.strip()
            fields = line.split("\t")
            tig = fields[0]
            if faidx[tig] < min_len:
                continue

            if tig not in breaks.keys():
                breaks[tig] = []
            # Reset at new tig
            if tig != prev_tig:
                if tig_start > 0:
                    breaks[prev_tig].append((tig_start, pos))

                prev_tig = tig
                split = False
                tig_start = 0
                break_start = 0

            pos, cov = int(fields[1]), float(fields[2])
            pos = pos - 1 # Because bed files start at coordinate 1, subtract 1 from pos

            # Initiate breaks only if coverage has been previously reached
            # in this contig, and we have reached beyond the end distance limit
            end_limit = faidx[tig] - distance
            if cov > min_cov and split == False and pos > distance and pos < end_limit:
                split = True

                # If there was a break_start position on this contig,
                # we have found a break
                if break_start > 0:
                    breaks[tig].append((tig_start, pos))
                    tig_start = break_start

            elif cov < min_cov and split:
                break_start = pos # Get start position for break
                split = False
    return { k:v for k,v in breaks.items() if v}

def readFasta(fasta):
    ''' Read fasta to dict.
    '''
    fasta_out = {} # sequences here
    faidx_out = {} # lengths here
    seq = []
    with open(fasta, "r") as fasta_in:
        for line in fasta_in:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    fasta_out[header] = "".join(seq)
                    faidx_out[header] = len(fasta_out[header])
                seq = []
                header = line.strip(">").split(" ")[0]

            else:
                seq.append(line)
        fasta_out[header] = "".join(seq)
        faidx_out[header] = len(fasta_out[header])
    return fasta_out, faidx_out

def writeBreaks(outprefix, breaks):
    with open(outprefix+".breaks.txt", "w") as breaksout:
        for k,v in breaks.items():
            for val in v:
                breaksout.write("{}\t{}\t{}\n".format(k,val[0],val[1]))

def breakFasta(fasta, breaks):
    fasta_out = {}
    for tig, seq in fasta.items():
        if tig not in breaks.keys():
            fasta_out[tig] = seq
        else:
            for idx, br in enumerate(breaks[tig]):
                start, end = br[0], br[1]
                fasta_out[tig+"."+str(idx)] = seq[start:end]

    return fasta_out

def writeFasta(outprefix, fasta):
    with open(outprefix+".fasta", "w") as out:
        for header, seq in fasta.items():
            out.write(">{}\n{}\n".format(header, seq))

def main():
    global faidx
    outprefix = args.out
    print("Reading fasta...")
    fasta, faidx = readFasta(args.fasta)
    print("Done.\nReading coverage...")

    breaks = readCov(args.bed, args.coverage, args.length, args.end_distance)
    print("Done.\nWriting breaks...")
    writeBreaks(outprefix, breaks)
    print("Done.\nBreaking fasta...")
    broken_fasta = breakFasta(fasta, breaks)
    print("Done.\nWriting fasta...")
    writeFasta(outprefix, broken_fasta)
    print("Done.")

if __name__ == "__main__":
    main()
