# -*- coding: utf-8 -*-
#python3 /media/fiji/yueming/20210410aQTL/bin/GTEx_V8_single_filter_SNPs_addcovariates.py -p /media/cuba/GTEx/temp/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP_2 -i /media/fiji/yueming/20210410aQTL/Input_Covariates/GTEx_covariates.v8.txt -o /media/fiji/yueming/20210410aQTL/filter/ >GTEx_V8_filter_SNPs.log 2>GTEx_V8_filter_SNPs.err &

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
def open_gz(filename):
	kind=filename.endswith(".gz")
	output = {}
	#ls=filename.split('/')
	if kind:
		gwas = gzip.open(filename,'rt', encoding='utf-8')
	else:
		gwas = open(filename)
	key = 1
	result = {}
	for line in gwas:
		if len(line) >0:
			if line[0:4] != "GTEX":			
				value = line.strip()
				if value == '0|0':
					value = '0'
				elif value == '0|1' or value == '1|0':
					value = '1'
				elif value == '1|1':
					value = '2'
				else:
					value = 'NA'
			else:
				value = line.strip()
				name = line.strip()
			result[key]=value
			key = key + 1
	gwas.close()
	return result,name
def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("-p","--pathfile",help="GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP_1")
	parser.add_argument("-o","--outputfile",help="output")
	parser.add_argument("-i","--inputfile",help="GTEx_covariates.v8.txt")
	args = parser.parse_args()
	pathfile = args.pathfile
	outputfile = args.outputfile
	inputfile = args.inputfile
	AllobFile = {}
	#covariates_dict = {}
	result_data,name = open_gz(pathfile)
	covariates = open(inputfile)
	for line in covariates:
		line=line.split('\t')
		if len(line) >2:
			GTEX = line[1].strip()
			#covariates_dict[GTEX] = line[1].strip()
			if GTEX == name:
				f_out = open(outputfile+line[0].strip()+"_"+name+".txt", "w")#创建文件对象
				for key1,valueL1 in result_data.items():
					f_out.write(''.join(valueL1) + '\n')
				f_out.close()		
if __name__ == '__main__':
	main(args)
