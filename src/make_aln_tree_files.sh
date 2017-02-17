#!/bin/bash

##This script will concatenate alignment fasta files and tree files for hyphy.
# if the sequences were simulated with codons and were translated to amino acids, make sure sim_space
# is set to 'codon'. This will allow to transform tree files into amino acid space by multiplying all 
# tree branches by 0.7702233
sim_space='codon' ##uncomment if the sequences were simulated in codon space and translated to amino acids
#sim_space='aa' ##uncomment if the sequences were simulated in amino acid space.

file_arr=(../substitution_matrices_in_pheno_models_data/translated_aln/*)

output_dir="hyphy/aln_tree_files/translated"

for aln_file in ${file_arr[*]} 
do	 
	temp1=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+'`
	tree_file="${temp1}.tre"
	temp2=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+_\d+'`
	aln_tree_file="${temp2}_aln_tree.txt"
	
	if [ $sim_space = 'codon' ]; then
		line=$(head -n 1 trees/$tree_file)
		t_codon_str=`echo $line | grep -oE '\d+\.\d+' | head -1`
		t_codon=bc <<< '$t_codon_str'
		#t_aa=$t_codon*0.7702233
		echo $t_codon
		#echo $t_aa
		#new_tree='(t1:'0.38',t2:'0.38');'
	fi
	#cat $aln_file trees/$tree_file > $output_dir/$aln_tree_file
done

