# -*- coding: utf-8 -*-
#python3 /media/disk3/CancerGWAS/CAUSALdb/run_fine_mapping.py -p /media/disk3/CancerGWAS/CAUSALdb/fine_mapping_input.txt >run_fine_mapping.log 2>run_fine_mapping.err &

#path	 -s	 -p

import argparse
import os
import sys
import gzip

args = sys.argv

def print_run(cmd):
	print(cmd)
	print("")
	os.system(cmd)
def AllFilePath(pathfilename):
	#root = os.getcwd() 
	file_names = open(pathfilename)
	file_ob_list = {}
	for file_name in file_names:
		file_name = file_name.split('\t')
		if file_name[0] != "path":
			fileob = file_name[0].strip()
			file_ob_list[fileob] = file_name[1].strip()+'\t'+file_name[2].strip()
	return file_ob_list
def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("-p","--pathfile",help="Fine_Mapping input pathfile")
	args = parser.parse_args()
	path = args.pathfile
	#AllPathFile = []
	output_dict = {}
	AllobFile = {}
	AllobFile = AllFilePath(path)
	for key,valueL in AllobFile.items():
		#valueL = ''.join(valueL)
		valueL = valueL.split('\t')
		filename = key.strip()
		if os.path.exists("./input/"+filename) and not os.path.exists("./output/"+filename.split(".")[0]+"_total_credible_set.txt") :#判断文件是否存在
			print_run("python3 fine_map_pipe.py -s " +valueL[0].strip()+ " -p " +valueL[1].strip()+ " ./input/"+filename+" output")
if __name__ == '__main__':
	main(args)


