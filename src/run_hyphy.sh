#!/bin/bash
num_sim=50
site_dupl=1
total_sites=130
taxa_num_arr=(128 256 512 1024 2048)
br_len_arr=(0.00005 0.0005 0.005 0.05 0.5)

if [ ! -d "./hyphy/rates" ]; then
	mkdir ./hyphy/rates
fi

if [ ! -d "./hyphy/rates/raw_rates" ]; then
	mkdir ./hyphy/rates/raw_rates
fi

model_name_arr=("JC" "JC_equalf" "WAG" "JTT" "LG")
for model in ${model_name_arr[@]}
do
	for br_len in `seq -f %.2f 0.02 0.02 0.52` 
	do	
		for n in $(seq 1 $num_sim) 
		do
			aln_tree_file=n2_bl${br_len}_${n}_aln_tree.txt
			rates_file=n2_bl${br_len}_${n}_${model}_rates.txt
			echo "INFILE = aln_tree_files/${aln_tree_file}" > hyphy/setup.txt
			echo "OUTFILE = rates/raw_rates/${rates_file}" >> hyphy/setup.txt
			echo "SITE_DUPL = ${site_dupl}" >> hyphy/setup.txt
			HYPHYMP hyphy/fitrates_${model}.bf
		done
	done
done
