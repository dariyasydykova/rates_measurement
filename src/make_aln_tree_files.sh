#!/bin/bash
file_arr=(./aln/site_dupl/*)

output_dir="hyphy/aln_tree_files/site_dupl"

for aln_file in ${file_arr[*]} 
do	 
	temp1=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+'`
	tree_file="${temp1}.tre"
	temp2=`echo $aln_file | grep -oE 'n\d+_bl\d+\.\d+_site_dupl_\d+_\d+'`
	aln_tree_file="${temp2}_aln_tree.txt"
	cat $aln_file trees/$tree_file > $output_dir/$aln_tree_file
done

