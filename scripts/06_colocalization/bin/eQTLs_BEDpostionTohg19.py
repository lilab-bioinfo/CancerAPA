# -*- coding: utf-8 -*-
#python /media/disk2/eQTL/bin/eQTLs_BEDpostionTohg19.py -p /media/disk2/eQTL/output/tempdata/bedQTLs_1 &

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

def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("-p","--patheqtl",help="eQTL pathfile")
	args = parser.parse_args()
	path_eqtl = args.patheqtl
	files_eqtl = os.listdir(path_eqtl)
	for file in files_eqtl:
		if not os.path.isdir(file):
			name = file.split('.')[0]
			rsidname = name+".txt"
			liftOverfile = open(path_eqtl+"/"+file)
			w = open(path_eqtl+"/../eQTLs_postionTohg19/"+rsidname,"w")
			w.write("chr"+"\t"+"position"+"\t"+"gene"+"\t"+"p-value"+"\t"+"N"+"\t"+"rsid"+"\n")
			for line in liftOverfile:
				#line = re.split('\t|_',line)
				line = line.split('\t')
				chromosome = line[0].strip()
				position = line[1].strip()
				other = line[3].strip().split('_')
				gene = other[0].strip()
				pvalue = other[1].strip()
				N = other[2].strip()
				rsid = other[3].strip()
				w.write(chromosome+"\t"+position+"\t"+gene+"\t"+pvalue+"\t"+N+"\t"+rsid+"\n")
			liftOverfile.close()
			w.close()
if __name__ == '__main__':
	main(args)
