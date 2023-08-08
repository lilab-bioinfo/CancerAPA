#zhec

#Use this script to get the .bed file which contains:
	#   V1   | V2  | V3  | V4               | V5        | V6
	#   chr# | TSS | TES | common_gene_name | ensemblID | strand

import sys
path_input=sys.argv[1]
path_output=sys.argv[2]
file_input=open(path_input,"r")
f_input = file_input.readlines()
file_input.close()
file_output=open(path_output,"w+")
i=0
gene=""
TSS=""
TES=""
while i<len(f_input):
	if i == 0 :	
		line_split=f_input[i].split('\t')
		chrNum=line_split[0]
		TSS=line_split[1]
		TES=line_split[2]
		gene=line_split[9].split(';')[4].split(' ')[2][2:-1]
		ensemblID=line_split[3]
		strand=line_split[5]
		#print(chrNum,TSS,TES,gene,ensemblID,strand)
		i=i+1
	else:
		line_split=f_input[i].split('\t')
		if gene == line_split[9].split(';')[4].split(' ')[2][2:-1]:
			if int(line_split[1]) < int(TSS):
				TSS=line_split[1]
			if int(line_split[2]) > int(TES):
				TES=line_split[2]
		else:
			file_output.write(chrNum+"\t"+TSS+"\t"+TES+"\t"+gene+"\t"+ensemblID+"\t"+strand+"\n")		
			chrNum=line_split[0]
			TSS=line_split[1]
			TES=line_split[2]
			gene=line_split[9].split(';')[4].split(' ')[2][2:-1]
			ensemblID=line_split[3]
			strand=line_split[5]
		i=i+1
		

file_output.close()
