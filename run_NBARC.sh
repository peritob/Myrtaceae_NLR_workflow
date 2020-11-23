#!/usr/bin/env bash

#to run this script: run_NBARC.sh interproscan_out.tsv augustus_out.gff3 > species_NBARC.gff3

# the first line of code produces a list to extract sequences from the braker amino acid fasta file
gawk 'BEGIN {FS="\t"} $6=="NB-ARC domain" { print $1 }' $1 | sort -k 1b,1 | uniq > NBARC_list

# this script outputs the gff3 file for all genes with NB-ARC domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.

gawk 'BEGIN {FS="\t"} $6=="NB-ARC domain" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp0

# then joins the gene id temp0 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all mathed lines and adds the scaffold location to column 4 and 5. The strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. This is likely to be the final true strandedness of the gene.  

join temp0 \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}'

rm temp0
