#!/bin/bash
num_sim=50
site_dupl=1
taxa_num_arr=(128 256 512 1024 2048)
br_len_arr=(0.00005 0.0005 0.005 0.05 0.5)
total_sites=129
sim_space='aa'

if [ -f "./src/run_sim_all_matrices.sh" ]; then
	rm ./src/run_sim_all_matrices.sh
fi

aln_dir=../substitution_matrices_in_pheno_models_data/aa_aln/all_sites
for br_len in ${br_len_arr[*]} 
do	 
	for num in ${taxa_num_arr[*]}  
	do
		for n in $(seq 1 $num_sim) 
		do
			tree_file=n${num}_bl${br_len}.tre
			aln_file=n${num}_bl${br_len}_${n}.fa
			echo python src/simulate_aln.py $sim_space trees/${tree_file} ${aln_dir}/${aln_file} ${site_dupl} ${total_sites} >> ./src/run_sim_all_matrices.sh	 
		done
	done
done

chmod +x ./src/run_sim_all_matrices.sh