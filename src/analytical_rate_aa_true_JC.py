import argparse
import textwrap
import sys
import os
import re
import numpy as np
from scipy import linalg

#the function derives Jukes-Cantor-like matrix with all substitution rates = 1/20 for all amino acids
def get_q_matrix():
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

#the function calculates a transition matrix P(r)
def get_p_matrix(t,r,q_matrix):
	p_matrix = linalg.expm(t*r*q_matrix) ##calculate P(t)=e^(t*r*Q)
	if p_matrix.sum()-20 > 0.0000001: #check that the rows add up to 1
		print("Rows in the projection matrix do not add up to 1!")
		sys.exit()
	return p_matrix

#the function generates rates pulled from a gamma distribution.
#Gamma distribution uses the shape parameter alpha=0.312
#and the rate parameter beta=1.027. 
#These were previously estimated for integrase HIV-1 proteins (Table 2 in Meyer and Wilke 2015)
def get_true_r_dict(m):
	true_r_dict=dict()
	true_r_lst=[]
	
	#pull one rate value from a gammadistribution
	for site in range(m):
		true_r_lst.append(np.random.gamma(0.312,1/1.027,1))
	
	#normalize rates
	for site in range(m):
		true_r_dict[site]=true_r_lst[site] / ( sum(true_r_lst)/m )
				
	return true_r_dict
	
def calculate_rate(outfile, m):
	
	#get true rates r for each site	
	true_r_dict=get_true_r_dict(m)
	
	q_matrix=get_q_matrix() #get the JC q matrix (not site specific)
	
	for t in np.arange(0.000002,2,0.02):
		r_dict=dict()
		for site in range(m):	
			true_r=true_r_dict[site]	
			p_matrix=get_p_matrix(t,true_r,q_matrix)
			
			r=0 #set up a counter to sum over p^(k)_i*p^(k)_ij(t) for all i,j 
			for i in range(20):
				for j in range(20):
	
					#skip summation when i = j
					if i==j:
						continue	
					#add when i does not equal j
					else:
						r+=(1.0/20.0)*p_matrix[i,j]
					
			r_dict[site]=np.log(1.0-(20.0/19.0)*r)
			
		#normalize site-wise rates by dividing them by the average rate in the sequence
		for site in r_dict:
			norm_r = r_dict[site] / ( sum(r_dict.values())/m )
			true_r = true_r_dict[site]
			
			line='%d,%f,%f,%f' %(site, t, true_r, norm_r) 
			outfile.write(line+'\n')

def main():

	'''
	Calculate analytical rate when the inference and true models uses the Jukes-Cantor-like matrix
	'''
	
	#creating a parser
	parser = argparse.ArgumentParser(
	formatter_class=argparse.RawDescriptionHelpFormatter,
			description='Calculate analytical rate when the inference and true models uses the Jukes-Cantor-like matrix',
	        epilog=textwrap.dedent('''\
	        Notation used in the description of columns r_tilde, r_tilde_small_t, and r_tilde_large_t:
	        r^(k)       - rate at site k
	        pi^(k)_i    - equilibrium frequency of amino acid i at site k
	        p^(k)_ij(t) - probability of an amino acid i changing into an amino acid j
	                      after time t at site k
	        m           - total number of sites
	        sum_ij      - sum over all amino acid i and sum over all amino acid j, 
	                      when i does not equal j
	        sum_l       - sum over all sites m
	        
            This script produces a CSV with the following columns: 
			
            Column name        Description
            =============================================================================
            site               Site in the amino acid sequence
            
            time               The time point at which the rate was calculated
            
            true_r             True rate (r) when the true model is P(r)=e^(t*r*Q)
            
            r_tilde            Site-wise rate defined by the equation 
                               r^(k) = log(1-(20/19)*sum_ij(pi^(k)_i*p^(k)_ij(t)) / 
                                       (1/m)*sum_l(log(1-(20/19)*sum_ij(pi^(l)_i*p^(l)_ij(t))))	
            '''))
	#adding arguments 
	parser.add_argument('m', metavar='<int>', type=int, help='number of sites to calculate rates')
	parser.add_argument('-o', metavar='<output.csv>', type=str, help='output rate file')
	args = parser.parse_args()
	
	if args.o is None:
		outfile='analytical_rates.csv'
	else:
		outfile=args.o
		
	#input arguments
	m=args.m
	
	##write a header for the output file
	outfile=open(outfile,'w')
	outfile.write('site,time,true_r,r_tilde\n')
	
	calculate_rate(outfile, m)

if __name__ == "__main__":
	main()     
