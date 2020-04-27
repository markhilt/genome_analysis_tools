#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Reads a fasta file and calculates N50.")
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

def calcN50(length_list):
    '''
    Calculate N50 from length_list
    '''
    length_list_sorted = sorted(length_list, reverse=True)
    half_length = sum(length_list) / 2
    n_contigs = 0
    tot_l = 0

    for l in length_list_sorted:
        n_contigs += 1
        tot_l += l

        if tot_l >= half_length:
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
    n50, l50 = calcN50(lengths)

    print("Total assembly size = {}\nNumber of scaffolds = {}\nN50 = {}\nL50 = {}\nLongest scaffold = {} bp".format(total, len(lengths), n50, l50, longest))

if __name__ == "__main__":
    main()
