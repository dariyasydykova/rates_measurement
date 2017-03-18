#!/bin/bash
site_dupl_arr=(10 100 1000 10000 100000)
num_sim=50
total_sites=129 #sets the limit. This will only generate dupl for one site
sim_space='aa'

if [ -f "./src/run_site_dupl_sim.sh" ]; then
	rm ./src/run_site_dupl_sim.sh 
fi

tree_file=n2_bl0.24.tre
aln_dir=../substitution_matrices_in_pheno_models_data/site_dupl_aln
for site_dupl in ${site_dupl_arr[*]}
do	
	for i in $(seq 1 $num_sim) 
	do
		aln_file=n2_bl0.24_site_dupl_${site_dupl}_${i}.fa
		echo python src/simulate_aln.py $sim_space trees/$tree_file ${aln_dir}/${aln_file} $site_dupl $total_sites >> ./src/run_site_dupl_sim.sh	
	done 
done

chmod +x ./src/run_site_dupl_sim.sh