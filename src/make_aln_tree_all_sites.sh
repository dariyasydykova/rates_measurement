#!/bin/bash
file_arr=($SCRATCH/substitution_matrices_in_pheno_models_data/aln/all_sites/*)

output_dir="hyphy/aln_tree_files/all_sites"

for aln_file in ${file_arr[*]} 
do	 
	temp1=`echo $aln_file | grep -oP 'n\d+_bl\d+\.\d+'`
	tree_file="${temp1}.tre"
	temp2=`echo $aln_file | grep -oP 'n\d+_bl\d+\.\d+_\d+'`
	aln_tree_file="${temp2}_aln_tree.txt"
	cat $aln_file trees/$tree_file > $output_dir/$aln_tree_file
done

