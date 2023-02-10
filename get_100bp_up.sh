#!/bin/bash
echo "extracting 100 bp upstream for ids in $1"
sed s'/\-RA//' $1 | 
grep -f - /projectnb/bradham/data/ReferenceSequences/wray-genome/L_var_clean.gff | 
awk '$3=="gene" { print }' | 
gff2bed | 
bedtools flank -i stdin -g ../output/sizes.genome -l 100 -r 0 -s |
bedtools getfasta -fi /projectnb/bradham/data/ReferenceSequences/wray-genome/Lvar_genome.fasta -bed stdin -name > 100bp_up.fa
