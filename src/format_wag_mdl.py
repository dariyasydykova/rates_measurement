import sys
import re
def reformat(tmpl, wag_mdl):
	for line in tmpl:
		line=line.strip()	
		
		if line.startswith("ModelMatrixName["):
			line_lst=line.split("=")
			idx = re.findall('\d+', line_lst[0])
			val = re.findall('\d*\.?\d+', line_lst[1])
		
			if "c" in line:
				if idx[0]=="14" and idx[1]=="0":
					wag_mdl.write("WAG_rt = {20,20,\n") 
				elif idx[0]=="19" and idx[1]=="17":
					wag_mdl.write("{"+idx[0]+","+idx[1]+",r*t*"+val[0]+"}};\n")			
				else:
					wag_mdl.write("{"+idx[0]+","+idx[1]+",r*t*"+val[0]+"}\n")
			else:
				if idx[0]=="19" and idx[1]=="17":
					wag_mdl.write("{"+idx[0]+","+idx[1]+",t*"+val[0]+"}\n")
				elif idx[0]=="14" and idx[1]=="0":
					wag_mdl.write("\n") 
					wag_mdl.write("WAG_t = {20,20,\n") 
				else:	
					wag_mdl.write("{"+idx[0]+","+idx[1]+",t*"+val[0]+"}\n")

		if line.startswith("equalFreqs["):
			line_lst=line.split("=")
			idx = re.findall('\d+', line_lst[0])
			val = re.findall('\d*\.?\d+', line_lst[1])
			if idx[0]=="0":
				wag_mdl.write("\n") 
				wag_mdl.write("WAG_freqs = {{"+val[0]+"}\n")	
			elif idx[0]=="17":
				wag_mdl.write("{"+val[0]+"}};")	
			else:
				wag_mdl.write("{"+val[0]+"}\n")			
def main():
	infile = sys.argv[1]
	outfile =  sys.argv[2]
	
	hyphy_wag_tmpl = open(infile,"r")
	wag_mdl = open(outfile,"w")
	
	reformat(hyphy_wag_tmpl,wag_mdl)

main ()