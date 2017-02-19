#!/bin/bash

##This script will concatenate alignment fasta files and tree files for hyphy.
file_arr=(../substitution_matrices_in_pheno_models_data/translated_aln/*)

output_dir="hyphy/aln_tree_files/translated"

for aln_file in ${file_arr[*]} 
do	 
	temp1=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+'`
	tree_file="${temp1}.tre"
	temp2=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+_\d+'`
	aln_tree_file="${temp2}_aln_tree.txt"
	cat $aln_file trees/$tree_file > $output_dir/$aln_tree_file
done

