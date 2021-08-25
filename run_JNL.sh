#!/usr/bin/env bash

#to run this script: run_JNL.sh interproscan_out.tsv augustus_out.gff3 > species_JNL.gff3

# to get a list of gene ids to make a fasta. 
#On the 5 August 2021 I removed $5=="G3DSA:1.10.8.430" from the LRR domain list as it is dubious.
#First line takes a tab delimited tsv file output from interproscan and searches column 5 for NBARC PfamID and prints the first column 'geneID' to file as unique.
#Second line does similar.
#Third line searches for any one of the column 5 IDs listed by OR function and prints the first column 'geneID' to file as unique. So any one or many of these may be present but if one- satisfies the search. Fourth line join adds tempU and tempW - if the geneID is present in both it prints out that ID to tempX.
gawk 'BEGIN {FS="\t"} $5=="PF00931" { print $1 }' $1 | sort -k 1b,1 | uniq > tempU
gawk 'BEGIN {FS="\t"} $5=="PF01419" { print $1 }' $1 | sort -k 1b,1 | uniq > tempV
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") { print $1 }' $1 | sort -k 1b,1 | uniq > tempW
join tempU tempV > tempX
join tempX tempW > JNL.list

# this script outputs the gff3 file for all genes with NB-ARC and Jacalin domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all Jacalin Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.

# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 

gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempU
gawk 'BEGIN {FS="\t"} $5=="PF01419" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempV
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempW
join tempU tempV > tempX
join tempX tempW > tempY

join tempY \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}'

rm tempU tempV tempW tempX tempY 
