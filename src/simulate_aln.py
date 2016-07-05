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
		site_num = int(line_lst[0])
		if site_num==num_model:
			break
		
		exp_temp = [ math.exp(-float(ddG)) for ddG in line_lst[1:]]
		pi_total = sum(exp_temp)
		pi_lst = [ exp/pi_total for exp in exp_temp]
		if 1-sum(pi_lst)>0.00001:
			print 'pi sum does not add up to 1'
		
		pi_dict[site_num] = pi_lst
	
	return pi_dict

def get_q_matrix(site_lst,num_model):
	'''
	function computes Q matrix for each site and returns the dictionary with Q matrix as the numpy array (values) for each site (keys).
	'''
	q_matrix_dict = {}
	
	for i in range(len(site_lst)):
		q_matrix_file = "q_matrices/site"+str(site_lst[i])+"_q_matrix_132L_A.txt"
		q_list = open(q_matrix_file,'r')
		q_matrix = []
		
		for line in q_list:
			line = line.strip()
			line_lst = line.split("\t")
			
			rates_list = [ float(x) for x in line_lst]
			q_matrix.append(rates_list)
	
		q_matrix_dict[site_lst[i]] = q_matrix
	
	return q_matrix_dict
			
def make_mutSel_model(site_lst,q_matrix_dict, pi_dict, tree_file, aln_file, site_dupl, num_model):
	'''
	function simulates MutSel model for the number of site_dupl with num_model model types
	'''
	tree=read_tree(file = tree_file)

	parts = []						
	for i in range(len(site_lst)):
		custom_matrix = q_matrix_dict[site_lst[i]]
		model = Model("custom", {"matrix": custom_matrix})
		
		
		##checking if frequencies are correct
		freqs=model.params['state_freqs']
		for j in range(len(freqs)):
			true_f = pi_dict[site_lst[i]]
			if freqs[j]-true_f[j]>0.01:
				print "wrong pi values"
		
		p = Partition(models = model, size = site_dupl)
		parts.append(p)
					
	##Evolving sequences	
	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = None, infofile = None, seqfile = aln_file)

def main(argv):

	if len(argv) != 6: # wrong number of arguments
		print '''Usage:
		simulate_aln.py <ddG_file> <tree_file> <aln_file> <site_dupl_num> <num_model>
		'''
		sys.exit()
		
	ddG_file = argv[1]
	tree_file = argv[2]
	aln_file = argv[3]
	
	# Define partition(s) 
	site_dupl = int(argv[4]) # number of sites simulated under one model
	num_model = int(argv[5]) 
	
	pi_dict = get_pi_dict(ddG_file,num_model)
	q_matrix_dict = get_q_matrix(pi_dict.keys(),num_model)
	make_mutSel_model(pi_dict.keys(),q_matrix_dict, pi_dict, tree_file, aln_file, site_dupl, num_model)
		
if __name__ == "__main__":
	main(sys.argv)