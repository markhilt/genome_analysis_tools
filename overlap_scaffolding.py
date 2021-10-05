#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""overlap_scaffoling.py

Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Description: Overlap-aware genome assembly reference-guided scaffolding. Takes
    a tab delimited file with contigs to merge, finds overlaps between them,
    and merges the contigs into a longer sequence. Borrows functions from
    ARBitR.

Copyright (c) 2021, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
"""

__version__ = "0.1"

import mappy as mp
import pysam
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", "--input_fasta", \
                    help="Input fasta file for contig merging.", \
                    type = str)
parser.add_argument("-p", "--paths", \
                    help="Input tsv file with contigs to merge.", \
                    type = str)
parser.add_argument("-o","--output", \
                    help="Prefix for output files.", \
                    type = str)
parser.add_argument("-a","--alignment_preset", \
                    help="Minimap2 alignment present. Valid values: 'map-ont', 'map-pb'", \
                    default = "map-ont", type = str)
parser.add_argument("-v","--version", help="Print version and quit.", action = "version", version = "overlap_scaffoling v.{}".format(__version__))
args = parser.parse_args()

def chop_cigar(cigar):
    '''
    Given a cigar string, chops it up into a list.
    '''
    cig_list = []
    pos = 0
    for idx, c in enumerate(cigar):
        if c.isalpha():
            cig_list.append(cigar[pos:idx+1])
            pos  = idx + 1
    return cig_list

def createConsensus(aln,string1,string2):
    '''Given two overlapping nucleotide strings and mappy alignment information,
    merges and returns the nucleotide strings. Overlap is assumed to be from suffix
    of string1 and prefix of string2

    Args:
        aln: <mappy.Alignment object>
        string1: first nucleotide string
        string2: second nucleotide string
    Returns:
        Str: The merged sequence.
    '''
    cig_list = chop_cigar(aln.cigar_str)
    # Track position in both strings
    ref_pos, query_pos = aln.r_st, aln.q_st
    output_string = [string1[:ref_pos]] # collect substrings in this list
    # Then walk through the aligned region in the cigar string,
    # gradually building the sequence. Because of the tendency of PacBio
    # to miss some bases, we will use the extra base at every indel position
    for cig in cig_list:
        if cig[-1] == "M":
            # If strings match, there is no problem
            # Use whatever sequence, they are the same
            # Then move both position trackers
            output_string.append( string1[ref_pos:ref_pos+int(cig[:-1])] )
            ref_pos += int(cig[:-1])
            query_pos += int(cig[:-1])
        elif cig[-1] == "I":
            # If insertion in query, add extra bases from query and increase position
            output_string.append( string2[query_pos:query_pos+int(cig[:-1])] )
            query_pos += int(cig[:-1])
        elif cig[-1] == "D":
            # If deletion in query, add extra bases from reference and increase position
            output_string.append( string1[ref_pos:ref_pos+int(cig[:-1])] )
            ref_pos += int(cig[:-1])

    # After iterating over the whole alignment, add remaining bases from
    # query sequence
    output_string.append( string2[query_pos:] )
    return ''.join(output_string)

def findOverlap(seq1,seq2, dist, alignment_preset):
    '''Find overlaps between two sequences.

    Description:
        Use mappy to find overlaps between all prefix-suffix pairs of the
        two given nucleotide sequences. The orientation is already determined
        so we consider only overlaps between the suffix of seq1 and prefix of
        seq2.

    Args:
        seq1 (str): First nucleotide sequence, serving as reference.
        seq2 (str): Second nucleotide sequence, serving as query.
        dist (int): Maximum distance from contig end to consider overlap
        alignment_preset (str): minimap2 preset

    Returns:
        list: suffix_overlaps, list of overlaps from the suffix of seq1
    '''
    #### First find reference suffix overlaps. We will treat seq1 as reference
    # Write seq1 to temporary fasta file
    pid = str(os.getpid())
    with open(pid+".tmp.fasta", "w") as tmpfasta:
        tmpfasta.write(">tmp\n")
        tmpfasta.write(seq1)

    # Build index from seq1
    #idx = mp.Aligner(fn_idx_in=pid+"tmp.fasta", preset="asm10")#k=19, w=10, scoring=[1,4,6,26])
    idx = mp.Aligner(fn_idx_in=pid+".tmp.fasta", preset=alignment_preset)

    # Align
    alignments = idx.map(seq2)

    # Iterate over alignments and search for overlapping ends
    suffix_overlaps = []
    if alignments:
        for aln in alignments:
            # Filter for alignments ending within the last dist bp of seq1,
            # and starting within the first 1 kb of seq2
            # We also need to control that the alignment is in the right direction
            # Overlap is at suffix of reference. Ideally it ends at the last
            # coordinate of reference, but because contig ends usually have
            # poor quality, the alignment might not reach that far. We can still
            # use the overlap, as long as there is an overhang from the query
            # sequence.
            if aln.r_en > aln.ctg_len - dist:
                # Overlap is at prefix of query, in this case it must be in
                # forward orientation
                if aln.q_st < dist and aln.strand == 1:
                    suffix_overlaps.append(aln)

    return suffix_overlaps

def merge_path(path):
    '''From a list of contigs, compute overlaps between each step and merge them.

    Args:
        path (list): list of contig names to be combined
            into a scaffold. The orientation of the contig should be appended
            at the end of the name as either f or r.

    Returns:
        superseq (str): Merged sequence.
    '''
    ctg = path[0][:-1] # Initiate at the first step
    dir = path[0][-1]
    superseq = fastafile[ctg] if dir == "f" else mp.revcomp(fastafile[ctg])
    # Then go through the list of steps in the path, align each step and merge
    for step in path[1:]:
        # Start by grabbing the contig name, orientation and seq
        ctg, dir = step[:-1], step[-1]
        seq = fastafile[ctg] if dir == "f" else mp.revcomp(fastafile[ctg])

        # Overlap with superseq.
        ovls = findOverlap(superseq,seq, 12000,args.alignment_preset)

        # If no overlap, merge by gap introduction and continue
        if not ovls:
            superseq += "N"*100
            superseq += seq
            continue
        # If more than one overlap, choose the one with the highest number of
        # aligned bases
        if len(ovls) > 1:
            ovl = ovls[0]
            for o in ovls:
                if ovl.mlen < o.mlen:
                    ovl = o
        else:
            ovl = ovls[0]
        # Next, merge the overlap
        superseq = createConsensus(ovl, superseq, seq)
    return superseq

def readPaths(pathsfile):
    with open(pathsfile, "r") as p:
        for line in p:
            line = line.strip()
            yield line.split("\t") # Generate lists of contigs in paths

def getOut():
    '''Creates a prefix for output files.
    '''
    if args.output:
        outfilename = args.output
    elif "/" in args.input_fasta:
        outfilename = args.input_fasta.split("/")[-1].split(".fasta")[0]+".merged"
    else:
        outfilename = args.input_fasta.split(".fasta")[0]+".merged"
    return outfilename

def writeSeqs(scaffolds):
    '''Write output
    '''
    with open(getOut()+".fasta", "w") as outfile:
        for idx, scf in enumerate(scaffolds, 1):
            outfile.write(">scaffold"+str(idx)+"\n")
            outfile.write(scf+"\n")

def main():
    global fastafile
    fastafile = pysam.FastaFile(args.input_fasta)
    scaffolds = []
    for path in readPaths(args.paths):
        # If only a single member in the path, the chromosome was already fully assembled
        # If so, only reverse comp if in reverse orientation.
        if len(path) == 1:
            ctg, dir = path[0][:-1], path[0][-1]
            seq = fastafile[ctg] if dir == "f" else mp.revcomp(fastafile[ctg])
            scaffolds.append(seq)
        else:
            # For each path, compute overlaps between member contigs,
            # and merge them into a longer sequence
            scaffolds.append(merge_path(path))
    fastafile.close()
    # Finally, write the new scaffolds to disk.
    writeSeqs(scaffolds)

if __name__ == "__main__":
    main()
