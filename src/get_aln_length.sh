#!/bin/bash
##this script checks simulated alignments lengths
aln_file=$1

cat $aln_file | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
