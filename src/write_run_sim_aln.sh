#!/bin/bash
num_sim=30
site_dupl=100000
total_sites=11
#sim_space='codon'
sim_space='aa'

if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

aln_dir=../substitution_matrices_in_pheno_models_data/aa_aln/ten_sites
for br_len in `seq -f %.2f 0.02 0.02 0.52` 
do	
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.fa
		
		echo python src/simulate_aln.py $sim_space trees/${tree_file} ${aln_dir}/${aln_file} ${site_dupl} ${total_sites} >> ./src/run_sim_aln.sh	 
	done
done

chmod +x ./src/run_sim_aln.sh