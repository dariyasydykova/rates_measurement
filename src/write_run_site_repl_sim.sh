#!/bin/bash
site_dupl_arr=(100 1000 10000 100000 1000000)
num_sim=50
total_sites=2 #sets the limit. This will only generate dupl for one site

if [ ! -d "./aln" ]; then
	mkdir aln
fi

if [ -f "./src/run_site_repl_sim.sh" ]; then
	rm ./src/run_site_repl_sim.sh 
fi

tree_file=n2_bl0.24.tre
for site_dupl in ${site_dupl_arr[*]}
do	
	for i in $(seq 1 $num_sim) 
	do
	aln_file=n2_bl0.24_site_dupl_${site_dupl}_${i}.fa
	echo python src/simulate_aln.py 132L_A_foldx_ddG.txt trees/$tree_file aln/site_dupl/$aln_file $site_dupl $total_sites >> ./src/run_site_repl_sim.sh	
	done 
done

chmod +x ./src/run_site_repl_sim.sh