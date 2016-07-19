#!/bin/bash
num_sim=50
site_dupl=1
total_sites=130
taxa_num_arr=(128 256 512 1024 2048)
br_len_arr=(0.00005 0.0005 0.005 0.05 0.5)


if [ ! -d "./hyphy/aln_tree_files/all_sites/" ]; then
	mkdir ./hyphy/aln_tree_files/all_sites/
fi

for br_len in ${br_len_arr[*]} 
do	 
	for num in ${taxa_num_arr[*]}  
	do
		for i in $(seq 1 $num_sim) 
		do
			tree_file=n${num}_bl${br_len}.tre
			aln_file=n${num}_bl${br_len}_${i}.fa
			aln_tree_file=n${num}_bl${br_len}_${i}_aln_tree.txt
			cat aln/all_sites/${aln_file} trees/${tree_file} > hyphy/aln_tree_files/all_sites/${aln_tree_file}
		done
	done
done

