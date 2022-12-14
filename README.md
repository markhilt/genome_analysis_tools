# genome_analysis_tools

This repo contains scripts that I frequenctly use in genome assembly analysis.

**find_telomeres.py** and **n50.py** are used for diagnosis of assemblies.

**asm_cov_break.py** and **convert_mol_cov.py** are used to parse output from tigmint molecule to break assemblies in regions where 10X barcode coverage drops.

**overlap_scaffoling.py** is used to build assembly scaffolds based on linkage information. Two files are needed: a fasta file with assembled sequences to build scaffolds from, and a tab delimited text file that contains the linkage information. This file should look like, e.g.:
```
contig1f  contig34r
contig2r  contig12f contig4f
contig5r
```
Each line represents a new scaffold. The contig names need to match the input fasta headers, except the final character that needs to be f or r for forward or reverse, to indicate the orientation that the contig should have in the scaffold. If a line contains a single contig name, only this contig will be put into the scaffold, but will be reverse complemented if the final character is r. overlap_scaffoling.py is overlap-aware, meaning that it will try to find overlaps between the ends of contigs before building the scaffold.

**merge_transcripts.py**: Use to merge sequences in a multiple sequence alignment in fasta format. Sequences are merged when they have the same header, up until "@". The header should contain e.g. a taxon name followed by @ and a sequence name. Script is useful when there are several representative sequences from a single taxon, e.g. because the gene is fragmented, and you want to have a single sequence per taxon for phylogenomics.
