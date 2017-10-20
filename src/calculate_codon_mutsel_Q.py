import argparse
import textwrap
import sys
import numpy as np

'''
This script calculates a codon (61x61) substitution matrix Q for the Mutation-Selection model. 
The elements of the matrix Q are q_ij defined as

q_ij=m_ij*(S_ij/(1-e^(-S_ij))), where m_ij is the mutation rate and S_ij=2*Ne*s_ij. 
Here, s_ij is the difference in fitness between codon i and codon j.
Ne is the effective population size. 

This script sets m_ij=1 for all codons i and j.
S_ij between codons i and j is set to equal S_ij between the amino acids that codon i and codon j belong to.

For example, the codon "TCT" codes for Serine and the codon for "CCT" codes for Proline. 
The S_ij between "TCT" and "CCT" was set to equal the S_ij between Serine and Proline.

Author: Dariya K. Sydykova
'''

##the codon list that will set the order for all codon matrices and vectors
ordered_codon_lst=['AAA','AAC','AAG','AAT',
	'ACA','ACC','ACG','ACT',
	'AGA','AGC','AGG','AGT',
	'ATA','ATC','ATG','ATT',
	'CAA','CAC','CAG','CAT',
	'CCA','CCC','CCG','CCT',
	'CGA','CGC','CGG','CGT',
	'CTA','CTC','CTG','CTT',
	'GAA','GAC','GAG','GAT',
	'GCA','GCC','GCG','GCT',
	'GGA','GGC','GGG','GGT',
	'GTA','GTC','GTG','GTT',
	'TAC','TAT', #stop codons TAA and TAG removed
	'TCA','TCC','TCG','TCT',
	'TGC','TGG','TGT', #stop codon TGA removed
	'TTA','TTC','TTG','TTT'
	]

#dictionary to help generate ddG values for the codons
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

###amino acid order in all amino acid matrices and vectors is specified by the ddG file
def get_ddg_dict(ddg_file):		

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
		
		##assign ddg values list to a site in a dictionary
		ddg_dict[site]=ddg_lst
	
	return ddg_dict, aa_lst #return the dictionary and the order of amino acids in the ddG file

#the function assigns amino acid ddG values to codons. 
def make_codon_ddg_dict(aa_ddg_dict, aa_lst):
		
	##The codon's ddG = ddG value of the amino acid that the codon belongs to 
	##codon ddG list will follow the codon order in ordered_codon_lst
	codon_ddg_dict = {} #set an empty dictionary to get codon ddG values for each site.
	for site in aa_ddg_dict:	
		aa_ddg_lst = aa_ddg_dict[site]
		codon_ddg_lst = np.empty(61)
		for i in range(len(aa_ddg_lst)):
			aa=aa_lst[i] #find the amino acid associated with the ddG values
			codon_subset = aa_to_codon_dict[aa] #find the codons that belonging to the amino acid
			for codon in codon_subset: #go through every codon 
				ind=ordered_codon_lst.index(codon) #find the index of the codon in the ordered list
				codon_ddg_lst[ind]=aa_ddg_lst[i] #insert ddG value to keep the codon order
		
		codon_ddg_dict[site]=codon_ddg_lst #assign the codon ddG values to a site	
	
	return codon_ddg_dict

def get_s_matrix(ddg_lst):

	s_matrix=np.empty((61,61)) #set up an empty matrix 61X61 to fill with S_ij values
	for i in range(61):
		for j in range(61):
			s_matrix[i,j]=ddg_lst[i]-ddg_lst[j] #fill the matrix at position i,j with S_ij=ddG_i-ddG_j
	
	return s_matrix

def get_q_matrix(s_matrix):
	##set q_ij=1 where the codon i and codon j are synonymous.
	q_matrix=np.empty((61,61))
	for i in range(61):
		for j in range(61):
			codon_i = ordered_codon_lst[i] 
			codon_j = ordered_codon_lst[j]
			
			##check if codon_i and codon_j belong to the same amino acid.
			for aa in aa_to_codon_dict:
				codon_subset=aa_to_codon_dict[aa]
				if codon_i in codon_subset and codon_j in codon_subset:
					syn=True
					break
				else:
					syn=False
			
			##assign 1 where two codons are synonymous and if they non-synonymous q_ij=S_ij/(1-e^(-S_ij))
			s_ij=s_matrix[i,j]
			if syn==True or s_ij==0:
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
	
	aa_ddg_dict, aa_lst =get_ddg_dict(ddg_file)
	codon_ddg_dict = make_codon_ddg_dict(aa_ddg_dict, aa_lst)
	for site in codon_ddg_dict:
		ddg_lst=codon_ddg_dict[site]
		s_matrix=get_s_matrix(ddg_lst)
		q_matrix=get_q_matrix(s_matrix)

		outfile='../q_matrix/codon/'+pdb_id+'_site%d_q_matrix.npy' %site
		np.save(outfile,q_matrix)
	
if __name__ == "__main__":
	main()     
