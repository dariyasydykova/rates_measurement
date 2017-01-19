import numpy as np
from scipy import linalg
import sys

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
			
		if site > site_limit:
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

##get_jc_q_matrix output a Jukes-Cantor like matrix with all the values not on the diagonal as 1/19 and values on the diagonal as -1. 
def get_jc_q_matrix():
	q_matrix = np.full((20,20),1.0/19)
	np.fill_diagonal(q_matrix,-1)
	return q_matrix
	
##get_p_matrix calculates a P(t) matrix for a given time with a given Q matrix.
##Here P(t) = e^(tQ) 
def get_p_matrix(t,q_matrix,r=1):
	p_matrix = linalg.expm(r*t*q_matrix)
	if np.sum(p_matrix)-20 > 0.0000001:
		print "Rows in the projection matrix do not add up to 1!"
		sys.exit()
		
	return p_matrix

##This functions makes a dictionary of true rates for each site provided. 
##The true rates are pulled from a gamma distribution. 
##The final rates are normalized to the mean rate amont all sites.
def get_true_r_dict(site_lst):
	true_r_dict = {}
	norm_r_dict = {}
	
	for site in site_lst:
		true_r = np.random.gamma(1,1)
		true_r_dict[site]=true_r
	
	mean_true_r = np.mean(true_r_dict.values())
	for site in site_lst:
		norm_r = true_r_dict[site]/mean_true_r
		norm_r_dict[site] = norm_r
	
	return norm_r_dict
	
##get_r_tilde calculates a site-wise rate (normalized) using equation ?? from D. K. Sydykova and C. O. Wilke (2017)
def get_r_tilde(site,t,ddg_dict,true_r_dict=False):
		
	##Calculate all sites denominator
	denom_sum = 0	
	for temp_site in ddg_dict:
		if true_r_dict==False:
			ddg_lst = ddg_dict[temp_site]
			pi_lst = get_pi_lst(ddg_lst)
			q_matrix = get_mutsel_q_matrix(ddg_lst)
			p_matrix = get_p_matrix(t,q_matrix)
		else:
			pi_lst = np.tile(1/20.0,20)
			q_matrix = get_jc_q_matrix()
			true_r = true_r_dict[temp_site]
			p_matrix = get_p_matrix(t,q_matrix,true_r)
		
		site_sum = 0
		for i in range(20):
			site_sum += pi_lst[i]*p_matrix[i,i]
			
		denom_sum += np.log((20/19.0)*site_sum-(1/19.0))
		
	##Calculate site-wise variables and the numerator
	if true_r_dict==False:
		site_ddg_lst = ddg_dict[site]
 		site_pi_lst = get_pi_lst(site_ddg_lst)
 		site_q_matrix = get_mutsel_q_matrix(site_ddg_lst)
  		site_p_matrix = get_p_matrix(t,site_q_matrix)
	else:
		site_pi_lst = np.tile(1/20.0,20)
		site_q_matrix = get_jc_q_matrix()
		site_true_r = true_r_dict[site]
		site_p_matrix  = get_p_matrix(t,q_matrix,site_true_r)
	
	site_sum = 0	
  	for i in range(20):
 		site_sum += site_pi_lst[i]*site_p_matrix [i,i]

 	#m is the total number of sites
 	m = len(ddg_dict.keys())
 	r_tilde = np.log( (20/19.0)*site_sum-(1/19.0) ) / ( (1.0/m) * denom_sum)
 	
 	if true_r_dict==False:
 		return r_tilde
 	else:
 		return true_r_dict,r_tilde
 
def get_r_tilde_small_t(site,ddg_dict):
	 
	##Calculate all sites denominator
	denom_sum = 0
	for temp_site in ddg_dict:		
		ddg_lst = ddg_dict[temp_site]
		pi_lst = get_pi_lst(ddg_lst)
		q_matrix = get_mutsel_q_matrix(ddg_lst)
		
		site_sum = 0
		for i in range(20):
			site_sum += pi_lst[i]*q_matrix[i,i]
		denom_sum += site_sum
		
	##Calculate site-wise variables and the numerator
	site_ddg_lst = ddg_dict[site]
 	site_pi_lst = get_pi_lst(site_ddg_lst)
 	site_q_matrix = get_mutsel_q_matrix(site_ddg_lst)
  	site_sum = 0
  	for i in range(20):
 		site_sum += site_pi_lst[i]*site_q_matrix [i,i]
	
 	#m is the total number of sites
 	m = len(ddg_dict.keys())
 	r_tilde_small_t = site_sum/ ( (1.0/m) * denom_sum)
 	
 	return r_tilde_small_t
 	
def main():

	if len(sys.argv) != 4: # wrong number of arguments
		print """Usage:
	python analytically_derived_rates.py <ddG_file> <output_txt_file> <site_limit> 
	"""
		sys.exit()

	infile = sys.argv[1]
	outfile = sys.argv[2]
	site_limit = sys.argv[3]
	
	if site_limit == "all":
		site_limit == None
	else:
		site_limit = int(site_limit)
	
	ddg_file = open(infile,"r")
	out_rate_file = open(outfile,"w")
	out_rate_file.write("site\ttime\tr_tilde_ms\tr_tilde_ms_small_t\tr_tilde_jc\ttrue_r_jc\n")
	
	ddg_dict = get_ddg_dict(ddg_file,site_limit)
	true_r_dict = get_true_r_dict(ddg_dict)
	
	for site in ddg_dict:				
		r_tilde_small_t = get_r_tilde_small_t(site,ddg_dict)
		for t in np.arange(0.000002,2,0.02):
			r_tilde_ms = get_r_tilde(site,t,ddg_dict)
			true_r_dict, r_tilde_jc = get_r_tilde(site,t,ddg_dict,true_r_dict) 
			
 			line = "%d\t%f\t%.10f\t%.10f\t%0.10f\t%0.10f\n" %(site,t,r_tilde_ms,r_tilde_small_t,r_tilde_jc,true_r_dict[site]) 
 			out_rate_file.write(line)
		
main()