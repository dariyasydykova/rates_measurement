#!/bin/bash

##This script will concatenate alignment fasta files and tree files for hyphy.

output_dir=trees/converted_trees

for br_len in `seq -f %.2f 0.02 0.02 0.52` 
do	
	tree_file=n2_bl${br_len}.tre
	python src/convert_codon_tree_to_aa_tree.py trees/$tree_file $output_dir
done