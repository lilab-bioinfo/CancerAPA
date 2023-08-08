import sys
import os

file1=sys.argv[1]
file2=sys.argv[2]
file3=sys.argv[3]
file1_open=open(file1,"r")
file2_open=open(file2,"r")
file3_write=open(file3,"w+")
file3_write.write("chr_pos chr.x position.x PV_eQTL N_eQTL rsid.x gene chr.y position.y rsid.y beta_GWAS varbeta_GWAS N_GWAS PV_GWAS MAF"+"\n")
dict_file1={}
f1_read=file1_open.readlines()
f2_read=file2_open.readlines()
gene=os.path.basename(file1)[:-4]
for i in range(1,len(f1_read)):
	line_split=f1_read[i].split('\t')
	key=line_split[0]+":"+line_split[1]
	value=f1_read[i]
	dict_file1[key]=line_split[0]+" "+line_split[1]+" "+line_split[2]+" "+line_split[3]+" "+line_split[4][:-1]
num=1
for i in range(1,len(f2_read)):
	line_split=f2_read[i].split('\t')
	chr_pos=line_split[2]
	if chr_pos in dict_file1:
		line=str(num)+" "+chr_pos+" "+dict_file1[chr_pos]+" "+gene+" "+line_split[0]+" "+line_split[1]+" "+line_split[3]+" "+line_split[4]+" "+line_split[5]+" "+line_split[6]+" "+line_split[7]+" "+line_split[8][:-1]
		file3_write.write(line+"\n")
		num +=1
	