#!/bin/bash

##This script will concatenate alignment fasta files and tree files for hyphy.
file_arr=(../../substitution_matrices_in_pheno_models_data/aa_aln/JC/*)
output_dir='../hyphy/aln_tree_files/JC'
#output_dir='hyphy/aln_tree_files/site_dupl'
#output_dir='hyphy/aln_tree_files/all_sites'

#sim_space='codon'
sim_space='aa'
site_dupl=false

for aln_file in ${file_arr[*]} 
do	 
	if [ $sim_space = 'codon' ]; then
		temp1=`echo $aln_file | grep -oE 'bl\d+\.\d+'`
		tree_file="converted_trees/n2_codon_${temp1}.tre"
	else
		temp1=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+'`
		tree_file="${temp1}.tre"
	fi	
	if [ $site_dupl = true ]; then
		temp2=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+_site_dupl_\d+_\d+'`
	else
		temp2=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+_\d+'`
	fi
	aln_tree_file="${temp2}_aln_tree.txt"
	echo $aln_tree_file
	cat $aln_file ../trees/$tree_file > $output_dir/$aln_tree_file
done

