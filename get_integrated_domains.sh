#!/usr/bin/env bash

prefix_name=$1
NLR_type=$2

#First extract all Pfam codes from the .tsv - thereby excluding all coils and other database outputs.

grep 'Pfam' ${prefix_name}_augustus_aa.fasta.tsv > ${prefix_name}_pfam.tsv

#Then extract all the genes based on the lists from NLR_clasification. For example,

grep -w -f ${prefix_name}_${NLR_type}.list ${prefix_name}_pfam.tsv > ${prefix_name}_${NLR_type}_Pfams.tsv

#Then extract all domains that are not present in NLR classifications using the parameter -v with fgrep. Make a list of NLR Pfams and call the txt file "NLR_Pfams.list" (PF00931, PF18052, PF08263, PF12799, PF13306, PF13855, PF13516, PF00560, PF07725, PF05659, PF01582, PF17862) associated with NLR genes.

grep -F -wv -f NLR_Pfams.list ${prefix_name}_${NLR_type}_Pfams.tsv > ${prefix_name}_${NLR_type}_ID_Pfams.tsv

#Obtain meaningful column data 

gawk 'BEGIN {FS="\t"; OFS="\t"}; {print $1, $5, $6, $7, $8}' ${prefix_name}_${NLR_type}_ID_Pfams.tsv > ${prefix_name}_${NLR_type}_ID_Pfams.data
