#!/bin/bash

# 2020 Benjamin J Perry and Steph Waller - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Development

printf "Starting getORF translation of filtered contigs...\n\n"
printf "Begin Execution at: $(date)\n\n"
sleep 1

###  Check for read files ###
if [ "$1" == "" ]; then
	printf "Error: No trimmed contig file passed\n"
  printf "Usage: getORF inputfile.fasta outputfile.faa\n\n"
	exit 1
fi

TODO: Fix the error message
if [ $(wc -l $1) -eq 0 ]; then
	printf "Error: No contigs in the file :( )\n"
	exit 1
fi

if [ "$2" == 0 ]; then
	printf "Error: No protein outfile named :( )\n"
  printf "Usage: getORF inputfile.fasta outputfile.faa\n\n"
	exit 1
fi

CONTIGS=$1
PROTEINS=$2

getorf -minsize 1200 -outseq $PROTEINS -sequence $CONTIGS

printf "Completed getorf translations.\n\n"
