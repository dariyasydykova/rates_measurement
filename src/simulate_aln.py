from pyvolve import *
import numpy as np
import sys
import random
import math
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def sim_mutSel_model(aln_type, tree_file, aln_file, site_dupl, site_limit):
	'''
	function simulates MutSel model for the number of site_dupl with num_model model types
	'''
	tree=read_tree(file = tree_file)

	parts = []	
	
	if aln_type=='aa':
		q_dir="../q_matrix/amino_acid/"
	elif aln_type=='codon':	
		q_dir="../q_matrix/codon/"
	else:
		print('wrong alignment type!')
		sys.exit()
	
	matrix_file_dict=dict() #store site list of 132L substitution matrix files
	q_matrices_files=os.listdir(q_dir)
	for q_matrix_file in q_matrices_files:
		match = re.search(r"site(\d+)_", q_matrix_file)
		site=int(match.group(1))
		matrix_file_dict[site]=q_matrix_file
	site_lst=sorted(matrix_file_dict.keys()) #extract a sorted list of sites
	final_site_lst=site_lst[:site_limit]
	
	for site in final_site_lst:
		q_matrix_file=matrix_file_dict[site]
		q_matrix=np.load(q_dir+q_matrix_file)

		model = Model("custom", {"matrix": q_matrix})		
		
		p = Partition(models = model, size = site_dupl)
		parts.append(p)
					
	##Evolving sequences	
	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = None, infofile = None, seqfile = aln_file)

def main(argv):

	if len(argv) != 6: # wrong number of arguments
		print('''Usage:
		simulate_aln.py <aa/codon> <tree_file> <aln_file> <site_dupl_num> <site_limit>
		''')
		sys.exit()
		
	aln_type = argv[1]
	tree_file = argv[2]
	aln_file = argv[3]
	
	# Define partition(s) 
	site_dupl = int(argv[4]) # number of sites simulated under one model
	site_limit = int(argv[5])
	
	sim_mutSel_model(aln_type, tree_file, aln_file, site_dupl, site_limit)
		
if __name__ == "__main__":
	main(sys.argv)