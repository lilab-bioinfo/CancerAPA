# -*- coding: utf-8 -*-
#module load R/3.6.2-anaconda3
#python /media/disk2/eQTL/bin/CrossMap_eQTLs_postionTohg19.py -p /media/disk2/eQTL/output/tempdata/eQTLs_rsid_new &
#python /media/disk2/eQTL/bin/CrossMap_eQTLs_postionTohg19.py -p /media/disk2/eQTL/output/tempdata/eQTLs_rsid_new_1 &
#python /media/disk2/eQTL/bin/CrossMap_eQTLs_postionTohg19.py -p /media/disk2/eQTL/output/tempdata/eQTLs_rsid_new_2 &

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
			gene = line[2].strip()
			pvalue = line[3].strip()
			N = line[4].strip()
			rsid = line[5].strip()
			temp_out.write(chromosome+"\t"+position+"\t"+position+"\t"+gene+"_"+pvalue+"_"+N+"_"+rsid+"\n")
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
	for file in files_eqtl:
		if not os.path.isdir(file):
			name = file.split('.')[0]
			bedname = name+".bed"
			temname = name+"_temp.txt"
			unlifted = bedname+".unmap"
			open_gz(path_eqtl+"/"+file,temname)
			#print_run("liftOver " + "tempfile.txt" + " " +  "/home/yueming/software/hg38ToHg19.over.chain.gz" + " " + "output.bed" + " " + unlifted)
			print_run("CrossMap.py bed ~/sQTL/input/hg38ToHg19.over.chain.gz "+ " " + temname + " "+ bedname)	
			print_run("rm " + " "  + temname)
			print_run("mv " + " " + unlifted+ " " + "unmapeQTLs")
			print_run("mv " + " " + bedname+ " " + "bedQTLs")
if __name__ == '__main__':
	main(args)
