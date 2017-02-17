import sys
import re

def main():

	if len(sys.argv) != 2: # wrong number of arguments
		print """Usage:
	python convert_codon_tree_to_aa_tree.py <codon_tree_file>
	"""
		sys.exit()

	tree_file = sys.argv[1]
	f = open(tree_file,'r')
	
	tree_top=f.readline()
	if re.search('', line):
		
	
main()
