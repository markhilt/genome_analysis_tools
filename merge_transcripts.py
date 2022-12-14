#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""merge_transcripts.py

Original author: Markus Hiltunen
GitHub: https://github.com/markhilt

Contributors: Markus Hiltunen

Description: Merge split sequences in an alignment.

Copyright (c) 2022
Licensed under the GPL3 license. See LICENSE file.
"""

__version__ = "0.1"

import argparse
import sys

parser = argparse.ArgumentParser( \
				description="	Merge split sequences in an alignment if the taxon \
						matches. Use when genes are split for some taxa. Taxon \
						and sequence names are delimited by @ in the fasta header.")
parser.add_argument("input", \
                    help="Input alignment multifasta.", \
                    type = str)
parser.add_argument("-v", "--version", help="Print version and quit.",
				action = "version",
				version = "mergeTranscripts v.{}".format(__version__))
args = parser.parse_args()

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def readFasta(multifastafile):
	''' Read input alignment file in fasta format.
	'''
	fa_dict = {}
	with open(multifastafile, "r") as multifasta:
		for line in multifasta:
			line = line.strip()
			if line.startswith(">"):
				taxon = line.strip(">")
				fa_dict[taxon] = ""
			else:
				fa_dict[taxon] = fa_dict[taxon] + line
	return fa_dict

def mergeSeqs(fasta):
	''' Merge sequences where taxa are the same.
	'''
	# Start by collecting all taxa where multiple sequences occur.
	taxa = [header.split("@")[0] for header in fasta.keys()]
	merged_dict = fasta
	def findDupl(l):
	    seen = set()
	    seen_add = seen.add
	    duplicates = set( x for x in l if x in seen or seen_add(x) )
	    return duplicates # Return a set of taxa that occur more than once
	for dupl in findDupl(taxa):
		dupl_headers = [ h for h in fasta.keys() if dupl in h]
		seq_len = len(fasta[dupl_headers[0]])
		merged_seq = list(fasta[dupl_headers[0]]) # Initial values from first sequence
		merged_dict.pop(dupl_headers[0])
		for h in dupl_headers[1:]: # Loop through remaining sequences
			seq = fasta[h]
			merged_dict.pop(h)
			for idx in range(seq_len):
				# If gap in new sequence, or position matches between seqs, do nothing
				if seq[idx] != "-" and seq[idx] != merged_seq[idx]:
					# If gap in old sequence, replace gap with new base/aa
					if merged_seq[idx] == "-":
						merged_seq[idx] = seq[idx]
					else:
						# Otherwise there is a conflict. Arbitrarily keep the
						# position from the initial sequence, but warn the user.
						eprint("WARNING: base conflict at sequence {} position {}!".format(h, str(idx)))
		merged_dict[dupl+"@merged"] = "".join(merged_seq)
	return merged_dict


def main():
	fa_dict = readFasta(args.input) # Read the alignment file into dict format.
	fa_merged = mergeSeqs(fa_dict)
	for k,v in fa_merged.items():
		print(">"+k)
		print(v)

if __name__ == '__main__':
	main()
