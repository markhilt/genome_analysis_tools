#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Reads a fasta file and calculates N50.")
parser.add_argument("-t", "--threshold", \
                    help="Threshold for output statistics (50 for N50 and L50) [50]", \
                    default = 50, \
                    type = int)
parser.add_argument("input_fasta", \
                    help="Input fasta file. Required.", \
                    type = str)
args = parser.parse_args()


def seqLengths(fastafile):
    '''
    Collect lengths of all sequences in fastafile
    '''
    with open(fastafile, "r") as fasin:
        lengths = []
        first = True
        seqlen = 0

        for line in fasin:
            line = line.strip()
            if line.startswith(">") and first == True:
                first = False

            elif line.startswith(">") and first == False:
                lengths.append(seqlen)
                seqlen = 0

            else:
                seqlen += len(line)
        lengths.append(seqlen)
    return lengths

def calcNX(length_list, threshold):
    '''
    Calculate NX and LX from length_list
    '''
    length_list_sorted = sorted(length_list, reverse=True)
    threshold_length = sum(length_list) * threshold * 0.01
    n_contigs = 0
    tot_l = 0

    for l in length_list_sorted:
        n_contigs += 1
        tot_l += l

        if tot_l >= threshold_length:
            return l, n_contigs

def formatGenomeSize(size):
    '''Reformat the size a human readable format
    '''
    if size >= 1000000000:
        # 1Gb
        return str(size / 1000000000) + " Gbp"
    elif size >= 1000000:
        # 1Mb
        return str(size / 1000000) + " Mbp"
    elif size >= 1000:
        # 1Mb
        return str(size / 1000) + " Kbp"

def main():
    lengths = seqLengths(args.input_fasta)
    longest = max(lengths)
    total = formatGenomeSize(sum(lengths))
    nX, lX = calcNX(lengths, args.threshold)

    print("Total assembly size = {0}\n \
    Number of scaffolds = {1}\n \
    N{2} = {3}\n \
    L{2} = {4}\n \
    Longest scaffold = {5} bp".format(total, len(lengths), args.threshold, nX, lX, longest))

if __name__ == "__main__":
    main()
