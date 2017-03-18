import numpy as np
from scipy import linalg
import sys

####################### Amino acid order in this scrip ####################
###amino acid order in all amino acid matrices and vectors: ACDEFGHIKLMNPQRSTVWY
def get_q_matrix():
	##set q matrix to JC
	q_matrix=np.empty((20,20))
	q_matrix.fill(1./20.)
		
	##set the rows to equal to zero by setting the diagonal to equal the negative of the sum of off-diagonal values
	np.fill_diagonal(q_matrix, 0) #set the diagonal to zero to be able to add rows together.
	np.fill_diagonal(q_matrix, -q_matrix.sum(axis=1))
	
	##check that each row adds up to zero
	if q_matrix.sum()-0 > 0.0000001:
		print("Rows in the subsitution matrix do not add up to 0!")
		sys.exit()
		
	return q_matrix

def get_true_r_dict(site_limit):
	true_r_dict={}
	true_r_lst=[]
	
	for i in range(site_limit):
		true_r=np.random.gamma(0.312,1/1.027,1)	
		true_r_lst.append(true_r)
		
	for i in range(site_limit):
		norm_r=true_r_lst[i]/np.mean(true_r_lst)
		true_r_dict[i]=norm_r
		
	return true_r_dict
	
def get_p_matrix(t,r,q_matrix):
	p_matrix = linalg.expm(r*t*q_matrix)
	if p_matrix.sum()-20 > 0.0000001:
		print("Rows in the projection matrix do not add up to 1!")
		sys.exit()
		
	return p_matrix
	
def get_pi_lst():
	pi_lst=np.empty(20)
	pi_lst.fill(1./20.)

	return pi_lst

def get_r_tilde(outfile, site_limit):
	
	##write a header for the output file
	outfile=open(outfile,'w')
	outfile.write('site\ttime\ttrue_r\tr_tilde\n')
	m=site_limit ##set total number of sites
	true_r_dict=get_true_r_dict(site_limit)
	
	for t in np.arange(0.000002,2,0.02):
		site_num_r=[]
		site_true_r=[]
		for site in range(site_limit):		
			q_matrix=get_q_matrix()
			true_r=true_r_dict[site]
			
			pi_lst = get_pi_lst()
			p_matrix=get_p_matrix(t,true_r,q_matrix)

			num_sum_r=0		
			for i in range(20):
				for j in range(20):
	
					##assign 0 where i = j
					if i==j:
						num_sum_r+=0
					else:
						num_sum_r+=pi_lst[i]*p_matrix[i,j]
					
			site_num_r.append(np.log(1.0-((20.0/19.0)*num_sum_r)))
			
		for i in range(len(site_num_r)):
			true_r = true_r_dict[i]
			site_r_tilde = site_num_r[i] / ( (1.0/m) * sum(site_num_r))

			line='%d\t%f\t%f\t%f' %(i+1, t, true_r, site_r_tilde) 
			outfile.write(line+'\n')

def main():

	if len(sys.argv) != 3: # wrong number of arguments
		print("""Usage:
		python calc_an_rates_aa_true_JC.py <output_txt_file> <site_limit> 
		""")
		sys.exit()

	outfile = sys.argv[1]
	site_limit = int(sys.argv[2])
	
	get_r_tilde(outfile, site_limit)

main()