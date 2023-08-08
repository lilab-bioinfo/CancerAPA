# -*- coding: utf-8 -*-
#cd ~/aQTL_pipeline
#python ~/aQTL_pipeline/bin/aQTL_format_PAINTOR_hg38.py -c ~/aQTL_pipeline/Whole_Blood.cis_eqtl_all_approachb.txt.gz -i ~/TWAS/input/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz -o ~/aQTL_pipeline/20210718_aQTLs_postionTohg38 >fastaqtl_gene_format.log 2>fastaqtl_gene_format.err &

#Author : Yueming Hu
#TIME:2021/7/07

import argparse
import os
import sys
import gzip
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr
#import scipy.stats
#from multiprocessing import Pool    # 导入进程池模块
import string
#使用导入的包函数
#stats = importr('stats')

args = sys.argv
gene_dict = {}
def print_run(cmd):
	print(cmd)
	print("")
	os.system(cmd)

def GetPath(pathfilename):
	#root = os.getcwd() #获取当前工作目录路径
	file_names = os.listdir(pathfilename)
	file_ob_list = {}
	for file_name in file_names:
		if len(file_name.split('.cis_eqtl_all_approachb.'))==2:
			fileob = pathfilename + '/' + file_name.strip() #循环地给这些文件名加上它前面的路径，以得到它的具体路径
		#file_ob_list[fileob] = file_name.split('.vcf.')[0]+ '.txt'
			file_ob_list[fileob] = file_name.strip()
	return file_ob_list
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
def open_gz(filename,outputfile,reference_dict,N):
	kind=filename.endswith(".gz")
	output = {}
	if kind:
		aQTL = gzip.open(filename,'rt', encoding='utf-8')
	else:
		aQTL = open(filename)
	first_line = aQTL.readline()  # 取第一行
	#first_line = first_line.split('\t')
	#key = 0
	output = {}
	#ls = filename.split('/')[-1].strip()
	f_out = open(outputfile, "w")#创建文件对象
	f_out.write("gene\tchr\tposition\trsid\tpval_nominal\tbeta\tvarbeta\tN\tMAF\n")
	for line in aQTL:
		line0 = line.split("\t")
		if len(line0)>=5:
			SNP = line0[0].strip()
			if SNP in reference_dict:
				gene = line0[1].strip()
				tstat = line0[3].strip()
				beta = line0[2].strip()
				if is_number(beta) and is_number(tstat):
					if tstat!="0":
						beta_se = str(float(beta)/float(tstat))
					else:
						beta_se = "0"
					valueL = reference_dict[SNP].split("\t")
					chrom = valueL[0].strip()
					position = valueL[1].strip()
					rsid = valueL[2].strip()
					MAF = valueL[3].strip()
					f_out.write(gene+"\t"+chrom+"\t"+position+"\t"+rsid+"\t"+line0[4].strip()+"\t"+beta+"\t"+beta_se+"\t"+str(N)+"\t"+MAF+"\n")
	aQTL.close()
	f_out.close()
	#return output
def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("-c","--aQTLfile",help="3'aQTL summarize_dap result Whole_Blood_gtex_v8.aqtl_annot.vcf.gz")
	parser.add_argument("-o","--outputfile",help="PAINTOR QTL input file")
	parser.add_argument("-i","--reference",help="1000G variant_metadata.txt.gz")
	#parser.add_argument("-t","--threads", type=int, default=1, help="number of threads [1]")
	args = parser.parse_args()
	pathfile = args.aQTLfile
	reference = args.reference
	outputfile = args.outputfile
	#threads = args.threads
	reference_dict = {}
	Allobfile = {}
	#Allobfile = GetPath(pathfile)
	# 创建一个最多开启threads进程的进程池
	kind=reference.endswith(".gz")
	if kind:
		variant = gzip.open(reference,'rt', encoding='utf-8')
	else:
		variant = open(reference)
	first_line = variant.readline()  # 取第一行
	for line in variant:
		line = line.split('\t')
		if len(line)>=7:
			chrom = "chr"+line[0].strip()
			position = line[1].strip()
			key = line[2].strip()
			rsid = line[6].strip()
			MAF = line[5].strip()
			if rsid =="NA":
				rsid = key
			reference_dict[key] = chrom+'\t'+position+'\t'+rsid+'\t'+MAF
	#po = Pool(threads)
	#for key1,valueL1 in Allobfile.items():
	name=pathfile.split("/")[-1].split('.cis_eqtl_all_approachb.')[0].strip()
	output=outputfile+"/"+name+ '.cis_aqtl_hg38.txt'
	if not os.path.exists(output): #判断文件是否存在
		temname = name+".cis_aqtl_hg19.txt"
		vcfaQTL = gzip.open("~/aQTL_pipeline/output_filter/"+name+"_SNPs.vcf_genotype.vcf.gz",'rt', encoding='utf-8')
		first_line = vcfaQTL.readline()  # 取第一行
		vcfaQTL.close()
		N = len(first_line.split('\t'))-1
		open_gz(pathfile,temname,reference_dict,N)
	# 		po.apply_async(open_gz,(key1,output,reference_dict))
	# print("----开始----")
	# # 关闭进程池,不再接收新的任务,开始执行任务
	# po.close()
	# # 主进程等待所有子进程结束
	# po.join()
	# print("----结束----")

if __name__ == '__main__':
	main(args)
