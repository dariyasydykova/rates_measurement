import sys
import re

def main():

	if len(sys.argv) != 3: # wrong number of arguments
		print """Usage:
	python convert_codon_tree_to_aa_tree.py <codon_tree_file> <output_dir>
	"""
		sys.exit()

	tree_file = sys.argv[1]
	output_dir = sys.argv[2]
	f = open(tree_file,'r')
	
	tree=f.readline()
	
	m = re.search('\d+\.\d+', tree)
	codon_bl=float(m.group(0))		
	aa_bl=codon_bl*0.7702233
	new_tree='(t1:'+str(aa_bl)+',t2:'+str(aa_bl)+');'
	
	new_tree_file=output_dir+'/'+'n2_codon_bl%.2f.tre' %codon_bl
	print new_tree_file
	out = open(new_tree_file,'w')
	out.write(new_tree)
	
main()
