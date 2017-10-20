import argparse
import textwrap
import sys
import numpy as np

'''
This script calculates an amino acid (20x20) substitution matrix Q for the Mutation-Selection model. 
The elements of the matrix Q are q_ij defined as

q_ij=m_ij*(S_ij/(1-e^(-S_ij))), where m_ij is the mutation rate and S_ij=2*Ne*s_ij. 
Here, s_ij is the difference in fitness between amino acid i and amino acid j
Ne is the effective population size. 

This script sets m_ij=1 for all amino acid i and j, and S_ij=ddG_i-ddG_j. 

Author: Dariya K. Sydykova
'''

###amino acid order in all amino acid matrices and vectors is: ACDEFGHIKLMNPQRSTVWY
def get_ddg_dict(ddg_file):		
	
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
		site = int(line_lst[0])	#record the site
		ddg_lst = line_lst[1:] 
			
		##rearrange ddG values to match the amino acid order in ordered_aa_lst
		ordered_ddg_lst=np.empty(20)
		for i in range(len(aa_lst)):
			aa = aa_lst[i]
			ind=ordered_aa_lst.index(aa)
			ordered_ddg_lst[ind]=ddg_lst[i]
		
		##assign rearranged ddg value lst to a site in a dictionary
		ddg_dict[site]=ordered_ddg_lst
	
	return ddg_dict

def get_s_matrix(ddg_lst):
		
	s_matrix=np.empty((20,20)) #set up an empty matrix 20X20 to fill with S_ij values
	for i in range(20):
		for j in range(20):
			s_matrix[i,j]=ddg_lst[i]-ddg_lst[j] #fill the matrix at position i,j with S_ij=ddG_i-ddG_j
	
	return s_matrix

def get_q_matrix(s_matrix):
	q_matrix=np.empty((20,20))
	for i in range(20):
		for j in range(20):
			
			##if i equals j, then q_ij=1
			##if i does not equal j, then q_ij=S_ij/(1-e^(-S_ij))
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

def main():
	'''
	Calculate amino acid substitution matrix Q using the mutation-selection model.
	'''
	
	#creating a parser
	parser = argparse.ArgumentParser(
	formatter_class=argparse.RawDescriptionHelpFormatter,
			description='Calculate amino acid substitution rates using ddG values file from Echave et al. Relationship between protein thermodynamic constraints and variation of evolutionary rates among sites. 2016',
	        epilog=textwrap.dedent('''\
            This script produces a .npy file containing a 2-D numpy array
            '''))
	#adding arguments 
	parser.add_argument('ddG', metavar='<ddG.txt>', type=str, help='input ddG file')
	parser.add_argument('pdb', metavar='<pdb_id>', type=str, help='PDB ID prefix')

	args = parser.parse_args()
	
	#input arguments
	pdb_id=args.pdb	
	ddg_file=args.ddG
	
	ddg_dict=get_ddg_dict(ddg_file)
	for site in ddg_dict:
		ddg_lst=ddg_dict[site]
		s_matrix=get_s_matrix(ddg_lst)
		q_matrix=get_q_matrix(s_matrix)
		
		outfile='../q_matrix/amino_acid/'+pdb_id+'_site%d_q_matrix.npy' %site
		np.save(outfile,q_matrix)
	
if __name__ == "__main__":
	main()     
