#!/bin/bash
# echo "$1"
echo "Grabbing ids from: $1"
seqtk subseq /projectnb/bradham/data/ReferenceSequences/wray-genome/L_var_transcripts.fasta $1 > regulons_transcripts.fa
