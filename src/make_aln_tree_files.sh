#!/bin/bash

##This script will concatenate alignment fasta files and tree files for hyphy.
file_arr=(../substitution_matrices_in_pheno_models_data/aa_aln/ten_sites/*)
output_dir="hyphy/aln_tree_files/ten_sites/"
#sim_space='codon'
sim_space='aa'

for aln_file in ${file_arr[*]} 
do	 
	temp1=`echo $aln_file | grep -oE 'bl\d+\.\d+'`
	if [ $sim_space = 'codon' ]; then
		tree_file="converted_trees/n2_codon_${temp1}.tre"
	else
		tree_file="n2_${temp1}.tre"
	fi	
	temp2=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+_\d+'`
	aln_tree_file="${temp2}_aln_tree.txt"
	cat $aln_file trees/$tree_file > $output_dir/$aln_tree_file
done

