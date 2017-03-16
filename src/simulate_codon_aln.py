from pyvolve import *
import numpy as np
import sys
import random
import math
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def sim_mutSel_model(tree_file, aln_file, site_dupl, site_limit):
	'''
	function simulates MutSel model for the number of site_dupl with num_model model types
	'''
	tree=read_tree(file = tree_file)

	parts = []	
	site_lst=[2,3,4,5,6,7,8,9,10,11]					
	for i in site_lst:
		q_matrix_file="q_matrices/codon/site%s_q_matrix_132L_A.txt" %i
		q_matrix=np.loadtxt(q_matrix_file)
		model = Model("custom", {"matrix": q_matrix})		
		
		p = Partition(models = model, size = site_dupl)
		parts.append(p)
					
	##Evolving sequences	
	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = None, infofile = None, seqfile = aln_file)

def main(argv):

	if len(argv) != 5: # wrong number of arguments
		print('''Usage:
		simulate_aln.py <tree_file> <aln_file> <site_dupl_num> <site_limit>
		''')
		sys.exit()
		
	tree_file = argv[1]
	aln_file = argv[2]
	
	# Define partition(s) 
	site_dupl = int(argv[3]) # number of sites simulated under one model
	site_limit = argv[4]
		
	if site_limit == "all" or site_limit=="None":
		site_limit == None
	else:
		site_limit = int(site_limit)
	
	sim_mutSel_model(tree_file, aln_file, site_dupl, site_limit)
		
if __name__ == "__main__":
	main(sys.argv)