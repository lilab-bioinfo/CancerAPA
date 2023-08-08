#!/bin/python3

#python3 /media/disk2/eQTL/bin/gtf_get_bed.py /media/disk2/eQTL/input/transcript/gencode.v31.annotation.gtf /media/disk2/eQTL/output/transcript/gencode.v31.annotation_gene.bed &

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
i=5
gene=""
TSS=""
TES=""
while i<len(f_input):
	if i == 5 :	
		line_split=f_input[i].split('\t')
		chrNum=line_split[0]
		TSS=line_split[3]
		TES=line_split[4]
		gene=line_split[8].split(';')[2].split(' ')[2][1:-1]
		ensemblID=line_split[8].split(';')[0].split(' ')[1][1:-1]
		strand=line_split[6]
		#print(chrNum,TSS,TES,gene,ensemblID,strand)
		i=i+1
	else:
		line_split=f_input[i].split('\t')
		if gene == line_split[8].split(';')[2].split(' ')[2][1:-1]:
			if int(line_split[3]) < int(TSS):
				TSS=line_split[3]
			if int(line_split[4]) > int(TES):
				TES=line_split[4]
		else:
			file_output.write(chrNum+"\t"+TSS+"\t"+TES+"\t"+gene+"\t"+ensemblID+"\t"+strand+"\n")		
			chrNum=line_split[0]
			TSS=line_split[3]
			TES=line_split[4]
			gene=line_split[8].split(';')[2].split(' ')[2][1:-1]
			ensemblID=line_split[8].split(';')[0].split(' ')[1][1:-1]
			strand=line_split[6]
		i=i+1
		

file_output.close()