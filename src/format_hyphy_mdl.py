import sys
import re
import numpy as np

def reformat(mat, freq, wag_mdl,model_name):
	
	wag_mdl.write(model_name+"_t = {20,20,\n")
	for i in range(20):
		for j in range(20):
			if i == j:
				continue	
			elif i==19 and j==18:
				wag_mdl.write("{"+str(i)+","+str(j)+",t*"+str(mat[i][j])+"}};\n")
				wag_mdl.write("\n")
				wag_mdl.write("\n")
			else:
				wag_mdl.write("{"+str(i)+","+str(j)+",t*"+str(mat[i][j])+"}\n")
				
	wag_mdl.write(model_name+"_rt = {20,20,\n")
	for i in range(20):
		for j in range(20):
			if i == j:
				continue	
			elif i==19 and j==18:
				wag_mdl.write("{"+str(i)+","+str(j)+",r*t*"+str(mat[i][j])+"}};\n")
				wag_mdl.write("\n")
				wag_mdl.write("\n")
			else:
				wag_mdl.write("{"+str(i)+","+str(j)+",r*t*"+str(mat[i][j])+"}\n")
	
	for i in range(20):
		if i==0:
			wag_mdl.write(model_name+"_freqs = {{"+str(freq[i])+"}\n")
		elif i==19:
			wag_mdl.write("{"+str(freq[i])+"}};")
		else:
			wag_mdl.write("{"+str(freq[i])+"}\n")
		

def collect_matrix_vals(tmpl):
	M = np.full((20,20),0)
	F = []
	
	j=0
	for line in tmpl:
		line=line.strip()
		line_lst=line.split(",")
			
		if "{" in line and "}" in line:
			if len(line_lst) == 20:
				val = re.findall('\d*\.?\d+', line)
		
				for i in range(len(val)):
					M[i,j] = float(val[i])
				j+=1
		
			if len(line_lst) == 1:
				val = re.findall('\d*\.?\d+', line)
				F=np.append(F,float(val[0]))
			
	return M,F
	
def main():
	infile = sys.argv[1]
	model_name = sys.argv[2]
	outfile =  sys.argv[3]
	
	hyphy_wag_tmpl = open(infile,"r")
	wag_mdl = open(outfile,"w")
	
	q_matrix, freq_arr = collect_matrix_vals(hyphy_wag_tmpl)
	reformat(q_matrix,freq_arr,wag_mdl,model_name)

main ()