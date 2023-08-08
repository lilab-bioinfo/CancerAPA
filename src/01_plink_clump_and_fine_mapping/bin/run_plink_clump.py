
# -*- coding: utf-8 -*-
#python $CUR_DIR/bin/run_plink_clump.py -p $CUR_DIR/bin/FinnGen_r5_numberfile.txt >run_plink_clump.log 2>run_plink_clump.err

# -s      -p      Class   CGID    Filename        File    GLGC
# 180989  EUR     Benign neoplasm CG0001  /media/Rome/yueming/20210125_GWASs/finngen_r5/finngen_R5_CD2_BENIGN_ADRENAL_EXALLC.gz   finngen_R5_CD2_BENIGN_ADRENAL_EXALLC.gz GLGC_CG0001_result.txt
# 180895  EUR     Benign neoplasm         CG0002  /media/Rome/yueming/20210125_GWASs/finngen_r5/finngen_R5_CD2_BENIGN_ANUS_ANAL_CANAL_EXALLC.gz   finngen_R5_CD2_BENIGN_ANUS_ANAL_CANAL_EXALLC.gz GLGC_CG0002_result.txt

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
		if file_name[0] != "-s":
			fileob = file_name[3].strip()+".txt"
			file_ob_list[fileob] = file_name[0].strip()+'\t'+file_name[1].strip()
	return file_ob_list
def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("-p","--pathfile",help="magma input pathfile")
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
		if valueL[1].strip() == "EUR":
			genome = "g1000_eur"
		elif valueL[1].strip() == "EAS":
			genome = "g1000_eas"
		elif valueL[1].strip() == "AMR":
			genome = "g1000_amr"
		elif valueL[1].strip() == "AFR":
			genome = "g1000_afr"
		else:
			genome = "g1000_eur"
		if os.path.exists("./plink_input/"+filename):
			if os.path.getsize("./plink_input/"+filename) > 55:
				if not os.path.exists("./plink_output/"+filename.split(".")[0]+".significant"):
					print_run("plink --clump ~/20210125_GWASs/output/plink_input/"+filename+" --bfile ~/software/magma_v1.08b_static/"+genome+"/"+genome+" --clump-r2 0.60 --clump-p1 0.00000005 --clump-p2 0.05 --out ~/20210125_GWASs/output/plink_output/"+filename.split('.')[0]+".significant")
				if not os.path.exists("./plink_output/"+filename.split(".")[0]+".lead"):
					print_run("plink --clump ~/20210125_GWASs/output/plink_input/"+filename+" --bfile ~/software/magma_v1.08b_static/"+genome+"/"+genome+" --clump-r2 0.10 --clump-p1 0.00000005 --clump-p2 0.05 --out ~/20210125_GWASs/output/plink_output/"+filename.split('.')[0]+".lead")
				if not os.path.exists("./plink_output/"+filename.split(".")[0]+".top"):
					print_run("plink --clump ~/20210125_GWASs/output/plink_input/"+filename+" --bfile ~/software/magma_v1.08b_static/"+genome+"/"+genome+" --clump-r2 0.10 --clump-p1 0.00000005 --clump-p2 0.05 --clump-kb 250 --out ~/20210125_GWASs/output/plink_output/"+filename.split('.')[0]+".top")				
if __name__ == '__main__':
	main(args)
