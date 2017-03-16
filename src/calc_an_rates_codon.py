import numpy as np
from scipy import linalg
import sys

####################### Amino acid and codon order in this scrip ####################
###amino acid order in all amino acid matrices and vectors: ACDEFGHIKLMNPQRSTVWY
###Codon order in all codon matrices and vectors:

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
	
	return ddg_dict,ordered_aa_lst
	
def get_codon_s_matrix(aa_ddg_lst,aa_lst,aa_to_codon_dict):
	
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
	
	##make a ddG codon list from amino acid ddG lst
	##codon ddG list will follow the codon order in ordered_codon_lst
	codon_ddg_lst=np.empty(61)
	for i in range(20):
		aa = aa_lst[i]
		codon_subset = aa_to_codon_dict[aa]
		for codon in codon_subset:
			ind=ordered_codon_lst.index(codon)
			codon_ddg_lst[ind]=aa_ddg_lst[i]
		
	s_matrix=np.empty((61,61))
	for i in range(61):
		for j in range(61):
			s_matrix[i,j]=codon_ddg_lst[i]-codon_ddg_lst[j]
	
	return s_matrix, ordered_codon_lst

def get_q_matrix(s_matrix, codon_lst, aa_to_codon_dict):
	##set q matrix[i,j]=1 where the codon i and codon j are synonymous.
	q_matrix=np.empty((61,61))
	for i in range(61):
		codon_i=codon_lst[i]
		for j in range(61):
			codon_j=codon_lst[j]
			
			
			##check if codon_i and codon_j belong to the same amino acid.
			for aa in aa_to_codon_dict:
				codon_subset=aa_to_codon_dict[aa]
				if codon_i in codon_subset and codon_j in codon_subset:
					syn=True
					break
				else:
					syn=False
			
			##assign 1 where two codons are synonymous and if they non-synonymous q_ij=s_ij/(1-e^(-s_ij))
			s_ij=s_matrix[i,j]
			if syn==True or s_ij==0:
				q_matrix[i,j]=1
			else:
				q_matrix[i,j]=s_ij/(1-np.exp(-s_ij))
	
	##set the rows to equal to zero by setting the diagonal to equal the negative of the sum of off-diagonal values
	np.fill_diagonal(q_matrix, -np.nansum(q_matrix,axis=1))
	
	##check that each row adds up to zero
	if q_matrix.sum()-0 > 0.0000001:
		print("Rows in the subsitution matrix do not add up to 0!")
		sys.exit()
		
	return q_matrix

def get_p_matrix(t,q_matrix):
	p_matrix = linalg.expm(t*q_matrix)
	if np.sum(p_matrix)-61 > 0.0000001:
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

def get_r_tilde(infile, outfile, site_limit, aa_to_codon_dict):
	
	##write a header for the output file
	outfile=open(outfile,'w')
	outfile.write('site\ttime\tr_tilde\tr_tilde_small_t\r_tilde_large_t\n')
	ddg_dict, aa_lst=get_ddg_dict(infile,site_limit)
	m=len(ddg_dict) ##set total number of sites
	
	for t in np.arange(0.000002,2,0.02):
		site_num_r=[]
		site_num_r_small_t=[]
		site_num_r_large_t=[]
		for site in ddg_dict:		
			ddg_lst = ddg_dict[site]
			s_matrix, codon_lst = get_codon_s_matrix(ddg_lst,aa_lst, aa_to_codon_dict)
			q_matrix=get_q_matrix(s_matrix, codon_lst, aa_to_codon_dict)

			q_matrix_file="q_matrices/codon/site%s_q_matrix_132L_A.txt" %site
			q_matrix.tofile(file=q_matrix_file,sep="\t")
	
			p_matrix_equil = get_p_matrix(10,q_matrix)
			pi_lst = get_pi_lst(p_matrix_equil)

			p_matrix=get_p_matrix(t,q_matrix)

			num_sum_r=0
			num_sum_r_small_t=0
			num_sum_r_large_t=0			
			for i in range(61):
				codon_i=codon_lst[i]
				for j in range(61):
					codon_j=codon_lst[j]
				
					##check if codon_i and codon_j belong to the same amino acid.
					for aa in aa_to_codon_dict:
						codon_subset=aa_to_codon_dict[aa]
						if codon_i in codon_subset and codon_j in codon_subset:
							syn=True
							break
						else:
							syn=False
	
					##assign 0 where two codons are synonymous 
					if syn==True:
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
			print(line)
			outfile.write(line+'\n')

def main():

	if len(sys.argv) != 4: # wrong number of arguments
		print("""Usage:
		python calc_an_rates_codon.py <ddG_file> <output_txt_file> <site_limit> 
		""")
		sys.exit()

	infile = sys.argv[1]
	outfile = sys.argv[2]
	site_limit = int(sys.argv[3])
	
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
	
	get_r_tilde(infile, outfile, site_limit, aa_to_codon_dict)

main()