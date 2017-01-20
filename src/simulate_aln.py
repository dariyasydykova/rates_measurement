from pyvolve import *
import numpy as np
import sys
import random
import math
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

##get_ddg_dict reformats delta-delta Gs from a ddg file. 
##The function outputs a dictionary that uses sites as keys and a list of delta delta Gs as values. 
def get_ddg_dict(ddg_file,site_limit):		
	ddg_file.next()
	ddg_dict = {}
	for line in ddg_file:
		line = line.strip()
		
		line_lst = np.fromstring(line,dtype=float,sep=' ')
		site = int(line_lst[0])	
		
		if site_limit == None:
			ddg_lst = np.delete(line_lst,0)
			ddg_dict[site]=ddg_lst
			continue
			
		if len(ddg_dict) == site_limit:
			break			
		else:
			ddg_lst = np.delete(line_lst,0)
			ddg_dict[site]=ddg_lst
	
	return ddg_dict

##get_pi_lst computes equilibrium frequencies (pi) for each amino acid from site-wise delta delta G values.
##The function uses equation 13 from J. Echave, E. J. Jackson, and C. O. Wilke (2015).
def get_pi_lst(ddg_lst):	
	pi_lst = np.exp(-ddg_lst)/np.sum(np.exp(-ddg_lst))
	if sum(pi_lst)-1 > 0.0000001:
		print "Equilibrium frequencies do not add up to 1!"
		sys.exit()
	return pi_lst
	
##get_q_matrix calculates a substitution matrix Q according to the Mutation-Selection model by Sella and Hirsh.
def get_mutsel_q_matrix(ddg_lst):	
	temp_s = np.reshape(ddg_lst, (len(ddg_lst), 1))
	s_matrix = temp_s - temp_s.transpose()
	np.fill_diagonal(s_matrix, 1)

	q_matrix=s_matrix/(1-np.exp(-s_matrix))
	np.fill_diagonal(q_matrix, 0)
	np.fill_diagonal(q_matrix, -np.nansum(q_matrix,axis=1))

	if q_matrix.sum()-0 > 0.0000001:
		print "Rows in the subsitution matrix do not add up to 0!"
		sys.exit()
		
	return q_matrix

def make_mutSel_model(ddg_file, tree_file, aln_file, site_dupl, site_limit):
	'''
	function simulates MutSel model for the number of site_dupl with num_model model types
	'''
	tree=read_tree(file = tree_file)

	ddg_dict = get_ddg_dict(ddg_file,site_limit)
	
	parts = []						
	for i in ddg_dict:
		ddg_lst = ddg_dict[i]
		custom_matrix = get_mutsel_q_matrix(ddg_lst)
		equilib_freqs = get_pi_lst(ddg_lst)
		model = Model("custom", {"matrix": custom_matrix})		
		
		##checking if frequencies are correct
		freqs=model.params['state_freqs']
		true_f = equilib_freqs
		for j in range(len(freqs)):
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
		simulate_aln.py <ddG_file> <tree_file> <aln_file> <site_dupl_num> <site_limit>
		'''
		sys.exit()
		
	ddg_file = argv[1]
	tree_file = argv[2]
	aln_file = argv[3]
	
	# Define partition(s) 
	site_dupl = int(argv[4]) # number of sites simulated under one model
	site_limit = int(argv[5]) 
	
	if site_limit == "all":
		site_limit == None
	else:
		site_limit = int(site_limit)
	
	ddg_file = open(ddg_file,"r")
	make_mutSel_model(ddg_file, tree_file, aln_file, site_dupl, site_limit)
		
if __name__ == "__main__":
	main(sys.argv)