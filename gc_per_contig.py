#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""gc_per_contig.py

Original author: Markus Thorén
GitHub: https://github.com/markhilt

Contributors: Markus Thorén

Description: Calculate the GC content per contig in a genome assembly.

Copyright (c) 2023
Licensed under the GPL3 license. See LICENSE file.
"""

__version__ = "0.1"

import argparse

parser = argparse.ArgumentParser( \
				description="Calculate GC content per contig")
parser.add_argument("input", \
                    help="Input multifasta.", \
                    type = str)
parser.add_argument("-v", "--version", help="Print version and quit.",
				action = "version",
				version = "gc_per_contig.py v.{}".format(__version__))
args = parser.parse_args()

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

def calcGC(fasta):
	''' Calculate GC content per contig.
	'''
	for header,seq in fasta.items():
		total_bases = len(seq)
		nGC = seq.count("G")+seq.count("C")+seq.count("c")+seq.count("g")
		GCpcnt = round(nGC/total_bases,3)
		print(header+"\t"+str(GCpcnt))

def main():
	fa_dict = readFasta(args.input) # Read the alignment file into dict format.
	calcGC(fa_dict)

if __name__ == '__main__':
	main()
