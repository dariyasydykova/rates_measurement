#!/bin/bash
aln_tree_file=n2_bl0.4_1_aln_tree.txt
rates_file=n2_bl0.4_1_rates.txt
echo "INFILE = aln_tree_input_files/${aln_tree_file}" > hyphy/setup.txt
echo "OUTFILE = rates/raw_rates/${rates_file}" >> hyphy/setup.txt
HYPHYMP hyphy/fitrates.bf
