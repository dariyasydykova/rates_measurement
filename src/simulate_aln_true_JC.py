from pyvolve import *
import numpy as np
import sys

def sim_aln(tree_file, aln_file, site_dupl):
    '''
    function simulates MutSel model for the number of site_dupl with num_model model types
    '''
    tree=read_tree(file = tree_file)

    parts = []  

    q_dir="../q_matrix/amino_acid/"

    matrix_file_dict=dict() #store site list of 132L substitution matrix files
    q_matrices_files=os.listdir(q_dir)
    print(q_matrices_files)
    for q_matrix_file in q_matrices_files:
        match = re.search(r"site(\d+)_JC_matrix", q_matrix_file)
        if match:
            site = int(match.group(1))
            matrix_file_dict[site] = q_matrix_file

    site_lst = sorted(matrix_file_dict.keys()) #extract a sorted list of sites

    for site in site_lst:
        q_matrix_file=matrix_file_dict[site]
        q_matrix=np.load(q_dir+q_matrix_file)

        model = Model("custom", {"matrix": q_matrix})       
        
        p = Partition(models = model, size = site_dupl)
        parts.append(p)
                    
    ##Evolving sequences    
    evolve = Evolver(partitions = parts, tree = tree)
    evolve(ratefile = None, infofile = None, seqfile = aln_file)

def main(argv):

    if len(argv) != 4: # wrong number of arguments
        print('''Usage:
        simulate_aln.py. <tree_file> <aln_file> <site_dupl_num>
        ''')
        sys.exit()
        
    tree_file = argv[1]
    aln_file = argv[2]
    
    # Define partition(s) 
    site_dupl = int(argv[3]) # number of sites simulated under one model
    
    sim_aln(tree_file, aln_file, site_dupl)
        
if __name__ == "__main__":
    main(sys.argv)