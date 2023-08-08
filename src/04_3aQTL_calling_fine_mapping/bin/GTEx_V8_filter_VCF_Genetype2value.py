# -*- coding: utf-8 -*-
#python3 ~/aQTL_pipeline/bin/GTEx_V8_filter_VCF_Genetype2value.py -p ~/aQTL_pipeline/output/Whole_Blood_SNPs.vcf_genotype.vcf -o ~/aQTL_pipeline/output_filter/Whole_Blood_SNPs.vcf_genotype.vcf >GTEx_V8_filter_SNPs.log 2>GTEx_V8_filter_SNPs.err &

#Author : Yueming Hu
#TIME:2021/04/12

import argparse
import os
import sys
import gzip
import pandas as pd
import numpy as np
args = sys.argv

def print_run(cmd):
	print(cmd)
	print("")
	os.system(cmd)

def GetPath(pathfilename):
	#root = os.getcwd() #获取当前工作目录路径
	file_names = os.listdir(pathfilename)
	file_ob_list = {}
	for file_name in file_names:  
		fileob = pathfilename + '/' + file_name #循环地给这些文件名加上它前面的路径，以得到它的具体路径
		#file_ob_list[fileob] = file_name.split('.')[0]
		file_ob_list[fileob] = file_name.strip()
	return file_ob_list
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
def transpose(matrix):
	new_matrix = []
	for i in range(len(matrix[0])):
		matrix1 = []
		for j in range(len(matrix)):
			matrix1.append(matrix[j][i])
		new_matrix.append(matrix1)
	return new_matrix
def open_gz(filename,outputfile):
	kind=filename.endswith(".gz")
	output = {}
	#ls=filename.split('/')
	if kind:
		gwas = gzip.open(filename,'rt', encoding='utf-8')
	else:
		gwas = open(filename)
	#result = {}
	f_out = open(outputfile, "w")#创建文件对象
	for line0 in gwas:
		line = line0.split('\t')
		if len(line) >2:
			key = line[0].strip()
			if key=="id":
				#result[key] = line0.strip()
				f_out.write(line0.strip() + '\n')
			else:
				k = 0
				for i in range(1,len(line)):
					if line[i]=='2':
						k+=1
				if k >=5:
					#result[key] = line0.strip()
					f_out.write(line0.strip() + '\n')
	gwas.close() 
	f_out.close() 
	#return result
def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("-p","--pathfile",help="GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz")
	parser.add_argument("-o","--outputfile",help="Input_Covariates/Covariates_by_tissue/*.v8.covariates.txt")
	#parser.add_argument("-o","--outputfile",help="output")
	args = parser.parse_args()
	pathfile = args.pathfile
	outputfile = args.outputfile
	#outputfile = args.outputfile
	AllobFile = {}
	#output_dict = {}
	result_data  = open_gz(pathfile,outputfile)
	#f_out = open(outputfile, "w")#创建文件对象
	#for key1,valueL1 in result_data.items():
	#	f_out.write(''.join(valueL1) + '\n')
	#f_out.close()		
if __name__ == '__main__':
	main(args)
