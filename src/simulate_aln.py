from pyvolve import *
import numpy as np
import sys
import random
import math
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_pi_dict(ddG_file,num_model):
	'''
	function computes pi vector for each site and returns a dictionary with pi vectors (values) for each site (keys). 
	'''
	ddG_list = open(ddG_file,'r')
	pi_dict = {}
	for line in ddG_list:
		if line.startswith('SITE'):
			continue
		
		line = line.strip()
		line_lst = line.split("\t")
		site_num = int(line_lst[0])-1
		if site_num==num_model+1:
			break
		
		exp_temp = [ math.exp(-float(ddG)) for ddG in line_lst[1:]]
		pi_total = sum(exp_temp)
		pi_lst = [ exp/pi_total for exp in exp_temp]
		if 1-sum(pi_lst)>0.00001:
			print 'pi sum does not add up to 1'
		
		pi_dict[site_num] = pi_lst
	
	return pi_dict
	
def make_mutSel_model(pi_dict, tree_file, aln_file,site_dupl,num_model):
	'''
	function simulates MutSel model for the number of site_dupl with num_model model types
	'''
	tree=read_tree(file = tree_file)

	parts = []						
	for i in range(num_model):
		aa_freq = pi_dict[i+1]
		model = Model("mutsel", {"state_freqs": aa_freq})
		p = Partition(models = model, size = site_dupl)
		parts.append(p)
					
	##Evolving sequences
	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = None, infofile = None, seqfile = aln_file)

def main(argv):

	if len(argv) != 4: # wrong number of arguments
		print '''Usage:
		simulate_aln.py <ddG_file><tree_file> <aln_file> 
		'''
		sys.exit()
		
	ddG_file = argv[1]
	tree_file = argv[2]
	aln_file=argv[3]

	# Define partition(s)
	site_dupl = 1000 # number of sites simulated under one model
	num_model = 5 #total number of models simulated such that the total number of sites is site_dupl*num_model
	
	pi_dict = get_pi_dict(ddG_file, num_model)
	make_mutSel_model(pi_dict,tree_file, aln_file,site_dupl,num_model)
		
if __name__ == "__main__":
	main(sys.argv)