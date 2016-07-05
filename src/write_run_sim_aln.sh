#!/bin/bash
num_sim=30
site_dupl=100000
total_sites=129

if [ ! -d "./aln" ]; then
	mkdir aln
fi

if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

for br_len in `seq -f %.2f 0.02 0.02 0.52` 
do	
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.fa
		echo "python src/simulate_aln.py 132L_A_foldx_ddG.txt trees/${tree_file} aln/${aln_file} ${site_dupl} ${total_sites}" >> ./src/run_sim_aln.sh	
	done
done

chmod +x ./src/run_sim_aln.sh