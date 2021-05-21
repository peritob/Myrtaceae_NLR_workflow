#!/usr/bin/env bash

#to run this script: run_NBLRR.sh interproscan_out.tsv augustus_out.gff3 > species_NBLRR.gff3

# to get a list of gene ids to make an NBSLRR.fasta
gawk 'BEGIN {FS="\t"} $5=="PF00931"  { print $1 }' $1 | sort -k 1b,1 | uniq > temp0
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:1.10.8.430" || $5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516" || $5=="PF12799" || $5=="PF13306" || $5=="PF13855")  { print $1 }' $1 | sort -k 1b,1 | uniq > temp1
join temp0 temp1 > NBLRR_list

# this script outputs the gff3 file for all genes with NB-ARC and LRR domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all LRR Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.

# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 

gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp0
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:1.10.8.430" || $5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp1
join temp0 temp1 > temp2

join temp2 \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}'

rm temp0 temp1 temp2
