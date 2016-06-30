import numpy as np
from scipy import linalg
import sys

def get_ddg_dict(ddg_file,total_sites):
	ddg_file.next()
	ddg_dict = {}
	for line in ddg_file:
		line = line.strip()
		
		line_arr = np.fromstring(line,dtype=float,sep=' ')
		site = int(line_arr[0])
		if site > total_sites:
			break
		else:	
			ddg_arr = np.delete(line_arr,0)
			ddg_dict[site]=ddg_arr
	return ddg_dict

def get_pi_dict(ddg_dict):
	pi_dict = {}
	for site in ddg_dict:	
		ddg_arr = ddg_dict[site]	
		pi_arr = np.exp(-ddg_arr)/np.sum(np.exp(-ddg_arr))
		pi_dict[site] = pi_arr
	return pi_dict
	
def get_q_matrix_dict_ms(ddg_dict):
	q_dict = {}
	for site in ddg_dict:
		ddg_arr = ddg_dict[site]
		temp_s = np.reshape(ddg_arr, (len(ddg_arr), 1))
		s_matrix = temp_s - temp_s.transpose()
		np.fill_diagonal(s_matrix, 1)
	
		q_matrix=s_matrix/(1-np.exp(-s_matrix))
		np.fill_diagonal(q_matrix, 0)
		np.fill_diagonal(q_matrix, -np.nansum(q_matrix,axis=1))
	
		if q_matrix.sum()-0 > 0.0000001:
			print "Rows in the subsitution matrix do not add up to 0!"
			break
		q_dict[site] = q_matrix
	
	return q_dict

def get_q_matrix_jc():
	q_matrix = np.full((20,20),1.0/19)
	np.fill_diagonal(q_matrix,-1)
	return q_matrix

def get_q_matrix_wag():
	infile = "./q_matrices/wag/wag.dat"
	q_file = open(infile,"r")
	
	q_matrix = []
	for line in q_file:
		line = line.strip()
		
		line_arr = np.fromstring(line,dtype=float,sep=' ')
		if line_arr[0]==-1:
			continue 
		
		if len(line_arr)==1:
			q_temp = np.hstack((line_arr,np.zeros(20-len(line_arr))))
		else:
			new_row = np.hstack((line_arr,np.zeros(20-len(line_arr))))
			q_temp=np.vstack((q_temp,new_row))
		
	q_tril_temp=np.matrix(q_temp)
	q_matrix=q_tril_temp+q_tril_temp.T
	np.fill_diagonal(q_matrix,q_tril_temp.diagonal())
	
	return q_matrix
		
def get_q_matrix_lg():
	q_matrix = np.full((20,20),1.0/19)
	np.fill_diagonal(q_matrix,-1)
	return q_matrix

def get_q_matrix_jtt():
	q_matrix = np.full((20,20),1.0/19)
	np.fill_diagonal(q_matrix,-1)
	return q_matrix
		
def get_p_matrix(t,q_matrix):
	p_matrix = linalg.expm(t*q_matrix)
	if np.sum(p_matrix)-20 > 0.0000001:
		print "Rows in the projection matrix do not add up to 1!"
	return p_matrix

def get_gamma_dist_true_r_lst(total_sites):
	true_r_lst = []
	for site in range(2,total_sites+1):
		true_r = np.random.gamma(1,1)
		true_r_lst.append(true_r)
	
	return true_r_lst
	
def get_site_sums(total_sites,pi_data,q_data,true_r_lst = None,site_specific = True):
	site_sum_dict = {}
	
	for site in range(2,total_sites+1):
		if site > total_sites:
			break
		else:  
			if site_specific:
				pi_arr = pi_data[site]
				q_matrix = q_data[site]
			else: 
				pi_arr = pi_data
				q_matrix = q_data
				
			true_r = true_r_lst[site-2]
			q_matrix = q_matrix*true_r
			
			sum_lst = []
			for t in np.arange(0.000002,2,0.02):
				p_matrix = get_p_matrix(t, q_matrix)
				sum = 0
				for i in range(20):
					prod = pi_arr[i]*p_matrix[i,i]
					sum += prod
				sum_lst.append(sum)
			site_sum_dict[site]= sum_lst
	
	return site_sum_dict
		
def get_r_tilde(total_sites,site_sum_dict):
	r_tilde_dict = {}
		
	sum_matrix=np.matrix(site_sum_dict.values())
	sum_denom_lst = sum_matrix.sum(axis=0)
	for site in site_sum_dict:
		r_tilde_lst = []
		for i in range(len(site_sum_dict[site])):
			denom_sum=sum_denom_lst[0,i]
			r_tilde = np.log((20.0/19.0)*site_sum_dict[site][i]-(1.0/19.0))/np.log((20.0/(19.0*(total_sites-1)))*denom_sum-(1.0/19.0))
			r_tilde_lst.append(r_tilde)
		r_tilde_dict[site] = r_tilde_lst
	return r_tilde_dict

def get_r_tilde_normalized(total_sites,site_sum_dict):
	r_tilde_dict = {}
		
	t_range = len(site_sum_dict[2])
	
	sum_lst = []
	for i in range(t_range):
		rate_sum = 0
		for site in site_sum_dict:		
			rate_sum += np.log((20.0/19.0)*site_sum_dict[site][i]-(1.0/19.0))	
		sum_lst.append(rate_sum)
	
	for site in site_sum_dict:
		r_tilde_lst = []
		for i in range(len(site_sum_dict[site])):
			r_tilde = ((total_sites-1)*np.log((20.0/19.0)*site_sum_dict[site][i]-(1.0/19.0)))/sum_lst[i]
			r_tilde_lst.append(r_tilde)
		r_tilde_dict[site] = r_tilde_lst
	return r_tilde_dict
	
def main():
	infile = sys.argv[1]
	outfile = sys.argv[2]

	ddg_file = open(infile,"r")
	out = open(outfile,"w")
	
	total_sites = 11
	ddg_dict = get_ddg_dict(ddg_file,total_sites)
		
	pi_dict_ms = get_pi_dict(ddg_dict)
	pi_dict_jc = np.full((20),1.0/20)
	q_dict_ms = get_q_matrix_dict_ms(ddg_dict)
	q_dict_jc = get_q_matrix_jc()
	q_dict_wag = get_q_matrix_wag() 
	#q_dict_lg = 
	#q_dict_jtt = 
	
	true_r_lst_ms = np.full(total_sites-1,1.0)
	true_r_lst_jc = get_gamma_dist_true_r_lst(total_sites)

	site_sums_ms = get_site_sums(total_sites, pi_dict_ms, q_dict_ms, true_r_lst_ms)
	site_sums_jc = get_site_sums(total_sites, pi_dict_jc, q_dict_jc, true_r_lst_jc, site_specific=False)
	site_sums_wag = get_site_sums(total_sites, pi_dict_jc, q_dict_wag, true_r_lst_jc, site_specific=False)

	r_tilde_dict_ms = get_r_tilde(total_sites,site_sums_ms)
	r_tilde_dict_jc = get_r_tilde(total_sites,site_sums_jc)
	r_tilde_dict_wag = get_r_tilde(total_sites,site_sums_wag)
	
	r_tilde_dict_ms_normalized =  get_r_tilde_normalized(total_sites,site_sums_ms)
	r_tilde_dict_jc_normalized = get_r_tilde_normalized(total_sites,site_sums_jc)
	r_tilde_dict_wag_normalized = get_r_tilde_normalized(total_sites,site_sums_wag)

	out.write("site\ttime\tr_tilde_ms_norm\tr_tilde_jc_norm\tr_tilde_wag_norm\ttrue_r_jc\n")
	t_lst = np.arange(0.000002,2,0.02)
	for site in ddg_dict:
		for i in range(len(t_lst)):
			line = str(site-1)+ "\t" +str(t_lst[i])+ "\t"+str(r_tilde_dict_ms_normalized[site][i])+"\t"+str(r_tilde_dict_jc_normalized[site][i])+"\t"+str(r_tilde_dict_wag_normalized[site][i])+"\t"+str(true_r_lst_jc[site-2])+"\n"
			out.write(line)
		
main()