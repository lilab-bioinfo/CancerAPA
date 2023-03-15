#!/bin/python
#Use this script to filter the eQTLs files and add the sample number of each tissue and change the gene name format.

import sys
path_eQTL=sys.argv[1]
samplenum=sys.argv[2]
path_output=sys.argv[3]

f1=open(path_eQTL,"r")
feQTL=f1.readlines()
f1.close()
f2=open(path_output,"w+")
f2.write("chr"+"\t"+"position"+"\t"+"gene"+"\t"+"p-value"+"\t"+"N"+"\n")
for line in feQTL:
	if len(line) == 0:
		continue
	line_split=line.split('\t')
	chrnum=line_split[1].split('_')[0]
	if chrnum == "variant":
		continue
	position=line_split[1].split('_')[1]
	gene=line_split[0]
	pvalue=line_split[6]
	#samplenum=line_split[3]
	f2.write(chrnum+"\t"+position+"\t"+gene+"\t"+pvalue+"\t"+samplenum+"\n")



