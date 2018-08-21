import argparse
import textwrap
import sys
import numpy as np

'''
This script calculates an amino acid (20x20) substitution matrix Q = r^(k)Q_JC, where r^(k) is the true rate of evolution at site k, and Q_JC is the Jukes-Cantor like matrix.

Q_JC has elements q_ij, defined as:

1/19, when i != j
1, when i = j

Author: Dariya K. Sydykova
'''
# the function extracts true site-wise rates used in analytical derivations
def get_true_rates(rate_file):
    file = open(rate_file, "r") # open the rate file to read
    rate_dict = dict() # create a dictionary to store rates for every site

    # read the file one line at a time
    for line in file:
        if line.startswith("site"): # skip the line with the header
            continue

        lst = line.split(",") # make a list of variables contained in the file
        site = int(lst[0]) # extract site
        rate = float(lst[2]) # extract rate

        # store site and its rate
        if site in rate_dict: # break the loop if sites start repeating
            break
        else: # add site-rate to the dictionary as a key-value pair
            rate_dict[site] = rate

    # return a dictionary with site-rate key value pairs
    return rate_dict

#the function derives Jukes-Cantor-like matrix with all substitution rates = 1/20 for all amino acids
def get_q_matrix():
	q_matrix=np.empty((20,20))
	q_matrix.fill(1./19.)
		
	##set the rows to equal to zero by setting the diagonal to equal the negative of the sum of off-diagonal values
	np.fill_diagonal(q_matrix, 0) #set the diagonal to zero to be able to add rows together.
	np.fill_diagonal(q_matrix, -q_matrix.sum(axis=0))
	
	##check that each row adds up to zero
	if q_matrix.sum()-0 > 0.0000001:
		print("Rows in the subsitution matrix do not add up to 0!")
		sys.exit()
		
	return q_matrix

def main():
    '''
    Calculate amino acid Jukes-Cantor-like substitution matrix.
    '''
    
    # creating a parser
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Calculate amino acid substitution rates using ddG values file from Echave et al. Relationship between protein thermodynamic constraints and variation of evolutionary rates among sites. 2016',
            epilog=textwrap.dedent('''\
            This script produces a .npy file containing a 2-D numpy array
            '''))
    # adding arguments 
    parser.add_argument('f', metavar='<rates.txt>', type=str, help='input file with true rates')

    args = parser.parse_args()
    
    # input arguments
    rate_file=args.f 

    # extract true rates for each site
    rate_dict = get_true_rates(rate_file)

    # loop over sites in the rates file
    for site in rate_dict:

        # create an amino acid Jukes-Cantor-like (JC-like) matrix
        q_JC = get_q_matrix()

        # extract true rate for a site
        rate = rate_dict[site]

        # calculate true substitution matrix, which is a product of true rate and JC-like matrix
        true_q = rate*q_JC

        # name output file that will store the matrix
        outfile='../q_matrix/amino_acid/site%d_JC_matrix.npy' %site

        # save the matrix to a .npy file
        np.save(outfile, true_q)
    
if __name__ == "__main__":
    main()     

