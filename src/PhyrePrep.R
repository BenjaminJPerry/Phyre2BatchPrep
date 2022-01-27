# bin/R
# 2022 Benjamin J Perry & Steph Waller - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1.0.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@otago.ac.nz
# Status: Development

# Description: This script takes and assembly from trinity and the virust summary out put table,
# it filters for contigs with no BLAST results and a minimum lenght and writes the corresponding
# contigs into a new file for translation prior to Phyre2 analysis.

# Usage: phyrePrep -in summary.txt -fast trinity.fasta -min 1000 -max ?

# Import libraries that we need
library(tidyverse) # dataframe manipution and IO
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html

# TODO : command line arg parsing

minSize = 1000
virusTableIn <- "sandbox/T1.summary.table.txt"
fastFileIn <- "sandbox/T1.trinity.fasta"
# TODO : implement minSize variable in the output file name.
outputFile <- "sandbox/filtered.trinity.len.1000.contigs.fasta"
# Read the virus summary table

# TODO: Remove command once argParse is setup!
setwd(dir = "~/Desktop/phyrePrep/")


virusTable <- read_tsv(file = virusTableIn, 
                       col_names = FALSE, 
                       trim_ws = TRUE, 
                       progress = TRUE,
                       col_types = "icidcccddccccddccdd")

# Read the transcriptome contigs (Trinity Assembly)


trinityContigs <- readDNAStringSet(filepath = fastFileIn,
                                   format = "fasta",
                                   use.names = TRUE)

# Filter the virus summary table for entries:
# 1. which have NO BLAST results
# 2. have length >= minSize
filteredVirusTable <- virusTable %>% filter(is.na(X7) & is.na(X13) & is.na(X17)) %>% filter(X3 >= minSize)

# extract contigs names from filteres virus summary table
filteredContigsNames <- filteredVirusTable$X2

# Take our filtered contigs names from virus summary table and
# select the corresponding fasta sequences

filteredContigs <- trinityContigs[filteredContigsNames]

# Print out the filtered contigs

writeXStringSet(x = filteredContigs,
                filepath = outputFile,
                append = FALSE,
                format = "fasta")

print("Thank you for using our awesome filtering tool XD K bye.")
