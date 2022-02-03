#!/bin/bash

# 2022 Benjamin J Perry and Steph Waller - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1.1
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@otago.ac.nz
# Status: Functional

printf "Starting getORF translation of filtered contigs...\n\n"
printf "Begin Execution at: $(date)\n\n"
sleep 1


###  Check files ###
if [ "$1" == "" ]; then
	printf "Error: No command line arguments given.\n\n"
  printf "Usage: getORF inputfile.fasta outputfile.faa LENGTH\n\n"
	exit 1
fi

# if [ "$1" == "" ]; then
# 	printf "Error: No trimmed contig file.\n\n"
#   printf "Usage: getORF inputfile.fasta outputfile.faa LENGTH\n\n"
# 	exit 1
# fi
# if [ "$2" == "" ]; then
# 	printf "Error: No protein output file given.\n\n"
#   printf "Usage: getORF inputfile.fasta outputfile.faa LENGTH\n\n"
# 	exit 1
# fi
# if [ "$3" == "" ]; then
# 	printf "Error: No minimum contig length given.\n\n"
#   printf "Usage: getORF inputfile.fasta outputfile.faa LENGTH\n\n"
# 	exit 1
# fi

CONTIGS=$1
PROTEINS=$2
LENGTH=$3

getorf -methionine No -maxsize ######## -minsize $LENGTH -outseq $PROTEINS -sequence $CONTIGS
