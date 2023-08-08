import sys
import os
path_input=sys.argv[1]
path_output=sys.argv[2]
f = open(path_input,"r")
file_input=f.readlines()
f.close()
i=0
gene=""
while i<len(file_input):
	if i == 0:
		line_split=file_input[1].split('\t')
		gene=line_split[2]
		w=open(path_output+"/"+gene+".txt","w+")
		title_split=file_input[i].split('\t')
		w.write(title_split[0]+"\t"+title_split[1]+"\t"+title_split[3]+"\t"+title_split[4]+"\t"+title_split[5])
		
	else:
		line_split=file_input[i].split('\t')
		if len(line_split) == 6:
			if gene == line_split[2]:
				w.write(line_split[0]+"\t"+line_split[1]+"\t"+line_split[3]+"\t"+line_split[4]+"\t"+line_split[5])	
			else:
				w.close()
				gene=line_split[2]
				w=open(path_output+"/"+gene+".txt","w+")
				title_split=file_input[0].split('\t')
				w.write(title_split[0]+"\t"+title_split[1]+"\t"+title_split[3]+"\t"+title_split[4]+"\t"+title_split[5])
				w.write(line_split[0]+"\t"+line_split[1]+"\t"+line_split[3]+"\t"+line_split[4]+"\t"+line_split[5])
	i=i+1
	if i == len(file_input):
		w.close()


