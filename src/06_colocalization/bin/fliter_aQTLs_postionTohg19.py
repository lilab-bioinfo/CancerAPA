# -*- coding: utf-8 -*-
#python ~/sQTL/bin/fliter_sQTLs_postionTohg19.py -p ~/sQTL/output/tempdata/sQTLs_rsid_new


#python CrossMap.py bed Whole_Blood_rsid.txt
#Author : Yueming Hu
#TIME:2020/12/03
'''
/media/disk2/eQTL/output/tempdata/eQTLs
'''
import argparse
import sys
import os
import re
import gzip
args = sys.argv
def print_run(cmd):
	print(cmd)
	print("")
	os.system(cmd)
def open_gz(filename,temname):
	kind=filename.endswith(".gz")
	#tempfile = {}
	#ls=filename.split('/')
	if kind:
		eQTL = gzip.open(filename,'rt', encoding='utf-8')
		#name = ls[5].split('.')[0]
	else:
		eQTL = open(filename)
		#name = ls[5]
	#first_line = eQTL.readline()  # 取第一行
	#first_line = first_line.split('\t')
	#key = 0
	temp_out = open(temname, "w")#创建文件对象
	for line in eQTL:
		line = line.split('\t')
		if line[0] != "chr":
			chromosome = line[0].strip()
			position = line[1].strip()
			#end = str(int(position)+1)
			intron = line[2].strip()
			pvalue = line[3].strip()
			N = line[4].strip()
			rsid = line[5].strip()
			temp_out.write(chromosome+"\t"+position+"\t"+position+"\t"+intron+"/"+pvalue+"/"+N+"/"+rsid+"\n")
			#tempfile[key] = chromosome+'\t'+position+'\t'+end+'\t'+chr_position+'\t'+rsid
			#key = key + 1
	eQTL.close()
	#return tempfile
def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("-p","--patheqtl",help="eQTL pathfile")
	args = parser.parse_args()
	path_eqtl = args.patheqtl
	files_eqtl = os.listdir(path_eqtl)
	print_run("mkdir -p " + " "+ "unmapeQTLs")
	print_run("mkdir -p " + " "+ "bedQTLs")
	print_run("mkdir -p " + " "+ "sQTLs_postionTohg19")
	for file in files_eqtl:
		if not os.path.isdir(file):
			name = file.split('.')[0]
			bedname = name+".bed"
			temname = name+"_temp.txt"
			unlifted = bedname+".unmap"
			open_gz(path_eqtl+"/"+file,temname)
			#print_run("liftOver " + "tempfile.txt" + " " +  "/home/yueming/software/hg38ToHg19.over.chain.gz" + " " + "output.bed" + " " + unlifted)
			print_run("CrossMap.py bed ~/sQTL/input/hg38ToHg19.over.chain.gz "+ " " + temname + " "+ bedname)			
			liftOverfile = open(bedname)
			w = open(path_eqtl+"/../sQTLs_intron_positionTohg19/"+file,"w")
			w.write("chr"+"\t"+"position"+"\t"+"intron"+"\t"+"p-value"+"\t"+"N"+"\t"+"rsid"+"\n")
			for line in liftOverfile:
				#line = re.split('\t|_',line)
				line = line.split('\t')
				chromosome = line[0].strip()
				position = line[1].strip()
				other = line[3].strip().split('/')
				intron = other[0].strip()
				pvalue = other[1].strip()
				N = other[2].strip()
				rsid = other[3].strip()
				w.write(chromosome+"\t"+position+"\t"+intron+"\t"+pvalue+"\t"+N+"\t"+rsid+"\n")
			liftOverfile.close()
			w.close()
			print_run("rm " + " "  + temname)
			print_run("mv " + " " + unlifted+ " " + "unmapeQTLs")
			print_run("mv " + " " + bedname+ " " + "bedQTLs")
if __name__ == '__main__':
	main(args)

