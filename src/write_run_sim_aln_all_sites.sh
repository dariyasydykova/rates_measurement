#!/bin/bash
num_sim=30
site_dupl=100000
total_sites=130
num_taxa_arr=(64 128 256 512)
br_len_arr=(0.00005 0.0005 0.005 0.05)

if [ ! -d "../substitution_matrices_in_pheno_models_data/aln" ]; then
	mkdir ../substitution_matrices_in_pheno_models_data/aln
fi

if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

for br_len in ${br_len_arr[*]} 
do	 
	for num in ${num_taxa_arr[*]}  
	do
		for i in $(seq 1 $num_sim) 
		do
			tree_file=n${num}_bl${br_len}.tre
			aln_file=n${num}_bl${br_len}_${i}.fa
			echo "python src/simulate_aln.py 132L_A_foldx_ddG.txt trees/${tree_file} ../substitution_matrices_in_pheno_models_data/aln/${aln_file} ${site_dupl} ${total_sites}" >> ./src/run_sim_aln.sh	
		done
	done
done

chmod +x ./src/run_sim_aln.sh