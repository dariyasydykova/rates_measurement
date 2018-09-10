import argparse
import textwrap
import sys
import os
import re
import numpy as np
from scipy import linalg

#the function calculates a transition matrix P(t)
def get_p_matrix(t,q_matrix): 
    p_matrix = linalg.expm(t*q_matrix) ##calculate P(t)=e^(t*Q)
    if p_matrix.sum()-20 > 0.0000001:#check that the rows add up to 1
        print("Rows in the projection matrix do not add up to 1!")
        sys.exit()
        
    return p_matrix
    
#the function calculates equilibrium frequencies from the transition matrix at equilibrium (t=10)
def get_pi_lst(q_matrix):
    p_matrix=get_p_matrix(10,q_matrix) #get the p matrix at the equilibrium
    pi_lst=np.diagonal(p_matrix) #the diagonal of the p matrix at the equilibrium is the equilibrium frequencies
    if np.sum(pi_lst)-1 > 0.01: #check that the equilibrium frequencies add up to 1
        print("Equilibrium frequencies do not add up to 1!")
        sys.exit()
    if np.all(pi_lst<=0): #check that the equilibrium frequences are not negative
        print('Negative equilibrium frequencies!')
        sys.exit()
    return pi_lst

#the function calculates site-wise rates
def calculate_rate(outfile, q_dir, m):
    
    #making sure the directory is specified correctly
    if q_dir.endswith("/"):
        pass
    else:
        q_dir+"/"   
    
    #get the list of sites to match m
    matrix_file_dict=dict() #store site list
    q_matrices_files=os.listdir(q_dir) #get all q matrices files in the directory q_dir
    for q_matrix_file in q_matrices_files:
        match = re.search(r"site(\d+)_", q_matrix_file)
        site=int(match.group(1))
        matrix_file_dict[site]=q_matrix_file
    site_lst=sorted(matrix_file_dict.keys()) #extract a sorted list of sites
    final_site_lst=site_lst[:m]
        
    for t in np.arange(0.000002,2,0.02):
        r_dict=dict()
        r_small_t_dict=dict()
        r_large_t_dict=dict()
        
        #loop over all site's p matrices and equilibrium frequencies to derive the numerator of the rate equations
        for site in final_site_lst:
            q_matrix_file=matrix_file_dict[site]
            q_matrix=np.load(q_dir+q_matrix_file)   #load the q matrix
            pi_lst = get_pi_lst(q_matrix) #calculate the equilibrium frequencies
            p_matrix=get_p_matrix(t,q_matrix) #calculate the transition matrix p

            r=0 #set up a counter to sum over p^(k)_i*p^(k)_ij(t) for all i,j 
            r_small_t=0 #set up a counter to sum over p^(k)_i*q^(k)_ij for all i,j 
            r_large_t=0 #set up a counter to sum over p^(k)_i*p^(k)_i for all i 
            for i in range(20):
                for j in range(20):
    
                    #skip summation when i = j
                    if i==j:
                        continue    
                    #add when i does not equal j
                    else:
                        r+=pi_lst[i]*p_matrix[i,j]
                        r_small_t+=pi_lst[i]*q_matrix[i,j]
                        r_large_t+=pi_lst[i]*pi_lst[j]
                    
            
            r_dict[site]=np.log(1.0-(20.0/19.0)*r)
            r_small_t_dict[site]=r_small_t
            r_large_t_dict[site]=np.log(1.0-(20.0/19.0)*r_large_t)
        
        #normalize site-wise rates by dividing them by the average rate in the sequence
        for site in r_dict:
            norm_r = r_dict[site] / ( sum(r_dict.values())/m )
            norm_r_small_t = r_small_t_dict[site] / ( sum(r_small_t_dict.values())/m )
            norm_r_large_t = r_large_t_dict[site] / ( sum(r_large_t_dict.values())/m )

            #write out the output rates
            line='%s,%f,%f,%f,%f' %(site, t, norm_r, norm_r_small_t, norm_r_large_t) 
            outfile.write(line+'\n')

def main():

    '''
    Calculate analytical rate when the true model is MutSel and the inference uses the Jukes-Cantor-like matrix
    '''
    
    #creating a parser
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Calculate analytical rate when the true model is MutSel and the inference uses the Jukes-Cantor-like matrix',
            epilog=textwrap.dedent('''\
            Notation used in the description of columns r_tilde, r_tilde_small_t, and r_tilde_large_t:
            r^(k)       - rate at site k
            pi^(k)_i    - equilibrium frequency of amino acid i at site k
            q^(k)_ij    - the substitution rate between amino acid i and amino acid j
            p^(k)_ij(t) - probability of an amino acid i changing into an amino acid j
                          after time t at site k
            m           - total number of sites
            sum_ij      - sum over all amino acid i and sum over all amino acid j, 
                          when i does not equal j
            sum_l       - sum over all sites m
            sum_i       - sum over all amino acid i
            
            This script produces a CSV with the following columns: 
            
            Column name        Description
            =============================================================================
            site               Site in the amino acid sequence
            
            time               The time point at which the rate was calculated
            
            r_tilde            Site-wise rate defined by the equation 
                               r^(k) = log(1-(20/19)*sum_ij(pi^(k)_i*p^(k)_ij(t)) / 
                                       (1/m)*sum_l(log(1-(20/19)*sum_ij(pi^(l)_i*p^(l)_ij(t))))
            
            r_tilde_small_t    Site-wise rate when t is small defined by the equation 
                               r^(k) = sum_ij(pi^(k)_i*q^(k)_ij) / 
                                       (1/m)*sum_ij(pi^(l)_i*q^(l)_ij)
            
            r_tilde_large_t    Site-wise rate when t is large defined by the equation 
                               r^(k) = log(1-(20/19)*sum_i((pi^(k)_i)^2) / 
                                       (1/m)*sum_l(log(1-(20/19)*sum_i((pi^(l)_i)^2)))      
            '''))
    #adding arguments 
    parser.add_argument('m', metavar='<int>', type=int, help='number of sites to calculate rates')
    parser.add_argument('-o', metavar='<output.csv>', type=str, help='output rate file')
    parser.add_argument('-q', metavar='<directory>', type=str, help='directory of the substitution matrix Q')

    args = parser.parse_args()
    
    if args.o is None:
        outfile='analytical_rates.csv'
    else:
        outfile=args.o
    
    if args.q is None:
        q_dir='../q_matrix/amino_acid/'
    else:
        q_dir=args.q
        
    #input arguments
    m=args.m
    
    ##write a header for the output file
    outfile=open(outfile,'w')
    outfile.write('site,time,r_tilde,r_tilde_small_t,r_tilde_large_t\n')
    
    calculate_rate(outfile, q_dir, m)

if __name__ == "__main__":
    main()     
