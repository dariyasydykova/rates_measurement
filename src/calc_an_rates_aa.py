import numpy as np
from scipy import linalg
import sys

####################### Amino acid order in this scrip ####################
###amino acid order in all amino acid matrices and vectors: ACDEFGHIKLMNPQRSTVWY

def get_ddg_dict(ddg_file,site_limit=None):		
	
	##the amino acid list that will set the order for all amino acid matrices and vectors
	ordered_aa_lst=['ALA','CYS','ASP','GLU', 'PHE',
	'GLY','HIS','ILE','LYS','LEU',
	'MET','ASN','PRO','GLN','ARG',
	'SER','THR','VAL','TRP','TYR']

	f = open(ddg_file,'r')
	ddg_dict = {} ##set an empty dictionary to get ddG values for each site.
	for line in f:
		line = line.strip()

		##get the order of amino acids in the ddg file
		if line.startswith('SITE'):
			aa_lst = line.split('\t')[1:]
			continue
			
		##get the site number and ddg values for each site 
		line_lst = np.fromstring(line,dtype=float,sep=' ')
		site = int(line_lst[0])	
		ddg_lst = line_lst[1:]
		
		##exit the loop if the file reached the final site
		if site > site_limit:
			break
			
		##rearrange ddg values to match the amino acid order in ordered_aa_lst
		ordered_ddg_lst=np.empty(20)
		for i in range(len(aa_lst)):
			aa = aa_lst[i]
			ind=ordered_aa_lst.index(aa)
			ordered_ddg_lst[ind]=ddg_lst[i]
		
		##assign rearranged ddg value lst to a site in a dictionary
		ddg_dict[site]=ordered_ddg_lst
	
	return ddg_dict
	
def get_codon_s_matrix(ddg_lst):
		
	s_matrix=np.empty((20,20))
	for i in range(20):
		for j in range(20):
			s_matrix[i,j]=ddg_lst[i]-ddg_lst[j]
	
	return s_matrix

def get_q_matrix(s_matrix):
	##set q matrix[i,j]=1 where the codon i and codon j are synonymous.
	q_matrix=np.empty((20,20))
	for i in range(20):
		for j in range(20):
			
			##assign 1 where i = j and if they non-synonymous q_ij=s_ij/(1-e^(-s_ij))
			s_ij=s_matrix[i,j]
			if i==j or s_ij==0:
				q_matrix[i,j]=1
			else:
				q_matrix[i,j]=s_ij/(1-np.exp(-s_ij))
	
	##set the rows to equal to zero by setting the diagonal to equal the negative of the sum of off-diagonal values
	np.fill_diagonal(q_matrix, 0) #set the diagonal to zero to be able to add rows together.
	np.fill_diagonal(q_matrix, -q_matrix.sum(axis=1))
	
	##check that each row adds up to zero
	if q_matrix.sum()-0 > 0.0000001:
		print("Rows in the subsitution matrix do not add up to 0!")
		sys.exit()
		
	return q_matrix

def get_p_matrix(t,q_matrix):
	p_matrix = linalg.expm(t*q_matrix)
	if p_matrix.sum()-20 > 0.0000001:
		print("Rows in the projection matrix do not add up to 1!")
		sys.exit()
		
	return p_matrix
	
def get_pi_lst(p_matrix):
	pi_lst=np.diagonal(p_matrix)
	if np.sum(pi_lst)-1 > 0.01:
		print("Equilibrium frequencies do not add up to 1!")
		sys.exit()
	if np.all(pi_lst<=0):
		print('Negative equilibrium frequencies!')
		sys.exit()
	return pi_lst

def get_r_tilde(infile, outfile, site_limit):
	
	##write a header for the output file
	outfile=open(outfile,'w')
	outfile.write('site\ttime\tr_tilde\tr_tilde_small_t\tr_tilde_large_t\n')
	ddg_dict=get_ddg_dict(infile,site_limit)
	m=len(ddg_dict) ##set total number of sites
	
	for t in np.arange(0.000002,2,0.02):
		site_num_r=[]
		site_num_r_small_t=[]
		site_num_r_large_t=[]
		for site in ddg_dict:		
			ddg_lst = ddg_dict[site]
			s_matrix = get_codon_s_matrix(ddg_lst)
			q_matrix=get_q_matrix(s_matrix)

			q_matrix_file="q_matrices/aa/site%s_q_matrix_132L_A.txt" %site
			np.savetxt(q_matrix_file, q_matrix)
	
			p_matrix_equil = get_p_matrix(10,q_matrix)
			pi_lst = get_pi_lst(p_matrix_equil)

			p_matrix=get_p_matrix(t,q_matrix)

			num_sum_r=0
			num_sum_r_small_t=0
			num_sum_r_large_t=0			
			for i in range(20):
				for j in range(20):
	
					##assign 0 where i = j
					if i==j:
						num_sum_r+=0
						num_sum_r_small_t+=0
						num_sum_r_large_t+=0	
					else:
						num_sum_r+=pi_lst[i]*p_matrix[i,j]
						num_sum_r_small_t+=pi_lst[i]*q_matrix[i,j]
						num_sum_r_large_t+=pi_lst[i]*pi_lst[j]
					
			site_num_r.append(np.log(1.0-((20.0/19.0)*num_sum_r)))
			site_num_r_small_t.append(num_sum_r_small_t)
			site_num_r_large_t.append(np.log(1.0-((20.0/19.0)*num_sum_r_large_t)))
		
		for i in range(len(site_num_r)):
			site_r_tilde = site_num_r[i] / ( (1.0/m) * sum(site_num_r))
			site_r_tilde_small_t = site_num_r_small_t[i] / ( (1.0/m) * sum(site_num_r_small_t))
			site_r_tilde_large_t = site_num_r_large_t[i] / ( (1.0/m) * sum(site_num_r_large_t))

			line='%d\t%f\t%f\t%f\t%f' %(i+1, t, site_r_tilde, site_r_tilde_small_t, site_r_tilde_large_t) 
			outfile.write(line+'\n')

def main():

	if len(sys.argv) != 4: # wrong number of arguments
		print("""Usage:
		python calc_an_rates_aa.py <ddG_file> <output_txt_file> <site_limit> 
		""")
		sys.exit()

	infile = sys.argv[1]
	outfile = sys.argv[2]
	site_limit = int(sys.argv[3])
	
	get_r_tilde(infile, outfile, site_limit)

main()