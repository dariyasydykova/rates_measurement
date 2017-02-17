import numpy as np
from numpy import linalg as LA
from scipy import linalg
import sys
import math

##get_ddg_dict reformats delta-delta Gs from a ddg file. 
##The function outputs a dictionary that uses sites as keys and a list of delta delta Gs as values. 
def get_ddg_dict(ddg_file,site_limit):		
	ddg_dict = {}
	for line in ddg_file:
		line = line.strip()
		
		if line.startswith('SITE'):
			aa_lst = line.split('\t')[1:]
			continue
			
		line_lst = np.fromstring(line,dtype=float,sep=' ')
		site = int(line_lst[0])	
		
		if site_limit == None:
			ddg_lst = np.delete(line_lst,0)
			ddg_dict[site]=ddg_lst
			continue
			
		if site > site_limit:
			break			
		else:
			ddg_lst = np.delete(line_lst,0)
			ddg_dict[site]=ddg_lst
	
	return ddg_dict,aa_lst

def get_nuc_mu_matrix(nuc,mu_rate):
	##Nucleotide mutation matrix
	
	nuc_lst = list(['T','G','C','A'])
	ind=nuc_lst.index(nuc)
	
	nuc_mu_matrix=np.ones((4,4))
	for i in range(4):
		for j in range(4):
			if i==ind:
				nuc_mu_matrix[i,j]=mu_rate
			elif j==ind:
				nuc_mu_matrix[i,j]=mu_rate
			else:
				continue	
	
	np.fill_diagonal(nuc_mu_matrix,0)
	np.fill_diagonal(nuc_mu_matrix, -np.nansum(nuc_mu_matrix,axis=1))
		
	if nuc_mu_matrix.sum()-0 > 0.0000001:
 		print "Rows in the nucleotide mutation matrix do not add up to 0!"
 		sys.exit()
 	
	return nuc_mu_matrix, nuc_lst
	
def get_aa_to_codon_dict():
	
	aa_to_codon_dict = {'PHE':['TTT','TTC'],
	'LEU':['TTA','TTG','CTT','CTC','CTA','CTG'],
	'ILE':['ATT','ATC','ATA'],
	'MET':['ATG'],
	'VAL':['GTT','GTC','GTA','GTG'],
	'SER':['TCT','TCC','TCA','TCG','AGC','AGT'],
	'PRO':['CCT','CCC','CCA','CCG'],
	'THR':['ACT','ACC','ACA','ACG'],
	'ALA':['GCT','GCC','GCA','GCG'],
	'TYR':['TAT','TAC'],
	'HIS':['CAT','CAC'],
	'GLN':['CAA','CAG'],
	'ASN':['AAT','AAC'],
	'LYS':['AAA','AAG'],
	'ASP':['GAT','GAC'],
	'GLU':['GAA','GAG'],
	'CYS':['TGT','TGC'],
	'TRP':['TGG'],
	'ARG':['CGT','CGC','CGA','CGG','AGA','AGG'],
	'GLY':['GGT','GGC','GGA','GGG']
	}
	return aa_to_codon_dict
	
def get_codon_s_matrix(ddg_lst,aa_lst):
	temp_s = np.reshape(ddg_lst, (len(ddg_lst), 1))
	aa_s_matrix = temp_s - temp_s.transpose()
		
	aa_to_codon_dict = get_aa_to_codon_dict()
	
 	codon_s_matrix = list()
 	codon_name_lst = list()
	for i in range(len(aa_s_matrix[0,:])):
		codon_s_row=list()
		
		for j in range(len(aa_s_matrix[:,0])):
			aa=aa_lst[j]
			codon_lst=aa_to_codon_dict[aa]
				
			aa_s_ij=aa_s_matrix[i,j]
			codon_s_ij_lst=[float(aa_s_ij)]*len(codon_lst)
			codon_s_row.extend(codon_s_ij_lst)
		
		aa=aa_lst[i]
		codon_lst=aa_to_codon_dict[aa]
		codon_name_lst.append(codon_lst)
		for k in range(len(codon_lst)):
			codon_s_matrix.append(codon_s_row)
	
	codon_s_matrix=np.array(codon_s_matrix)
 	return codon_s_matrix,codon_name_lst
	
##convert a nucleotide mutation rate matrix to a codon mutation rate matrix
def get_codon_mut_rate(nuc_mu_matrix, nuc_lst,codon_lst_lst):
	
	nuc_mu_dict={}
	for i in range(4):
		for j in range(4):
			nuc_mu=nuc_mu_matrix[i,j]
			ni_nj=nuc_lst[i]+nuc_lst[j]
			nuc_mu_dict[ni_nj]=nuc_mu
	
	codon_lst=[codon for sublst in codon_lst_lst for codon in sublst]
	
	codon_mu_matrix=np.zeros((61,61))
	for i in range(61):	
		for j in range(61):
			ci = codon_lst[i]
			cj = codon_lst[j]
			pos_diff = [k for k in range(len(ci)) if ci[k] != cj[k]]
			
			if len(pos_diff)==1:
				change_nuc = ci[pos_diff[0]]+cj[pos_diff[0]]
				codon_mu_rate=nuc_mu_dict[change_nuc]
				codon_mu_matrix[i,j]=codon_mu_rate
			else:
				codon_mu_matrix[i,j]=0
	
	np.fill_diagonal(codon_mu_matrix, -np.nansum(codon_mu_matrix,axis=1))
			
	if np.sum(codon_mu_matrix)-61 > 0.0000001:
		print "Rows in the mutation rate matrix do not add up to 0!"
		sys.exit()
		
	return codon_mu_matrix

##get_q_matrix calculates a substitution matrix Q according to the Mutation-Selection model by Sella and Hirsh.
def get_mutsel_q_matrix(s_matrix,mu_matrix):	
		
	q_matrix=np.ones((61,61))
	non_zero=s_matrix != 0
	q_matrix[non_zero]=mu_matrix[non_zero]*(s_matrix[non_zero]/(1-np.exp(-s_matrix[non_zero]))) ##A_ij=m_ij*(S_ij/(1-exp(-S_ij))
	
	np.fill_diagonal(q_matrix,0)
	np.fill_diagonal(q_matrix, -np.nansum(q_matrix,axis=1))

 	if q_matrix.sum()-0 > 0.0000001:
 		print "Rows in the subsitution matrix do not add up to 0!"
 		sys.exit()

 	return q_matrix
	
##get_p_matrix calculates a P(t) matrix for a given time with a given Q matrix.
##Here P(t) = e^(tQ) 
def get_p_matrix(t,q_matrix):
	p_matrix = linalg.expm(t*q_matrix)
	if np.sum(p_matrix)-61 > 0.0000001:
		print "Rows in the projection matrix do not add up to 1!"
		sys.exit()
		
	return p_matrix
	
def get_pi_lst(p_matrix):
	pi_lst=np.diagonal(p_matrix)
	if np.sum(pi_lst)-1 > 0.01:
		print "Equilibrium frequencies do not add up to 1!"
		sys.exit()
	if np.all(pi_lst<=0):
		print 'Negative equilibrium frequencies!'
		sys.exit()
	return pi_lst

def get_br_len_factor(ddg_dict,aa_lst,nuc_mu_matrix,nuc_lst):

	f_num=0
	f_denom=0
	f_site_sum_num=0
	f_site_sum_denom=0
	for site in ddg_dict:
		ddg_lst = ddg_dict[site]
	
		s_matrix, codon_lst_lst = get_codon_s_matrix(ddg_lst,aa_lst)
		mu_matrix=get_codon_mut_rate(nuc_mu_matrix, nuc_lst, codon_lst_lst)
		q_matrix=get_mutsel_q_matrix(s_matrix,mu_matrix)
				
		p_matrix_equil = get_p_matrix(10,q_matrix)
		pi_lst = get_pi_lst(p_matrix_equil)
		
		N_matrix=np.ones((61,61))
		N_matrix[s_matrix==0]=0
		
		for i in range(len(pi_lst)):
			for j in range(len(pi_lst)):
				if i==j:
					continue
				else:
					mu_ij=mu_matrix[i,j]
					N_ij=N_matrix[i,j]
					pi=pi_lst[i]
					f_num += mu_ij*N_ij*pi
					f_denom += mu_ij*pi
		
		f_site_sum_num += f_num
		f_site_sum_denom += f_denom
	
	f=f_site_sum_num/f_site_sum_denom
	return f	

def main():

	if len(sys.argv) != 5: # wrong number of arguments
		print """Usage:
	python calculate_br_len_factor.py <ddG_file> <site_limit> <nuc> <mu_rate> 
	"""
		sys.exit()

	infile = sys.argv[1]
	site_limit = sys.argv[2]
	nuc_mu_rate = sys.argv[3]
	mu_rate = float(sys.argv[4])
	
	if site_limit == "all":
		site_limit == None
	else:
		site_limit = int(site_limit)
	
	ddg_file = open(infile,"r")
	
	ddg_dict,aa_lst = get_ddg_dict(ddg_file,site_limit)

	nuc_mu_matrix, nuc_lst = get_nuc_mu_matrix(nuc_mu_rate,mu_rate)
	f = get_br_len_factor(ddg_dict,aa_lst,nuc_mu_matrix,nuc_lst)
	print f
	
main()