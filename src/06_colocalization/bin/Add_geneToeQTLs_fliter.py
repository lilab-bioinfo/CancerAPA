# -*- coding: utf-8 -*-
#python /media/disk2/eQTL/bin/Add_geneToeQTLs_fliter.py -b /media/disk2/eQTL/output/transcript/gencode.v31.annotation_gene.bed -p /media/disk2/eQTL/output/tempdata/eQTLs &
#python /media/disk2/eQTL/bin/Add_geneToeQTLs_fliter.py -b /media/disk2/eQTL/output/transcript/gencode.v36.annotation_gene.bed -p /media/disk2/eQTL/output/tempdata/eQTLs &
#python /media/disk2/eQTL/bin/Add_geneToeQTLs_fliter.py -b /media/disk2/eQTL/output/transcript/hg38_stable_geneID_TO_name.txt -p /media/disk2/eQTL/output/tempdata/eQTLs_rsid &

#Author : Yueming Hu
#TIME:2020/12/03
'''
/media/disk2/eQTL/output/tempdata/eQTLs
'''
import argparse
import sys
import os

args = sys.argv
def main(args):
	parser = argparse.ArgumentParser()
	#parser.add_argument("-b","--bedfile",help="hg38 gencode.v31.annotation_gene.bed pathfile")
	parser.add_argument("-b","--bedfile",help="hg38_stable_geneID_TO_name.txt pathfile")
	parser.add_argument("-p","--patheqtl",help="eQTL pathfile")
	args = parser.parse_args()
	path_eqtl = args.patheqtl
	bedfile = args.bedfile
	bed_Dict = {}
	files_eqtl = os.listdir(path_eqtl)
	bed = open(bedfile)
	for line in bed:
		line = line.split("\t")
		if line[0].strip()!="Gene name":
			gene_id = line[1].strip()
			bed_Dict[gene_id] = line[0].strip()
	for file in files_eqtl:
		if not os.path.isdir(file):
			f = open(path_eqtl+"/"+file)
			w = open(path_eqtl+"/../sQTLs_rsid_new/"+file,"w")
			w.write("chr"+"\t"+"position"+"\t"+"gene"+"\t"+"p-value"+"\t"+"N"+"\t"+"rsid"+"\n")
			for line in f:
				line = line.split("\t")
				if line[0].strip()!="chr":
					CHROM = line[0].strip()
					position = line[1].strip()
					gene_id = line[2].strip().split(".")[0]
					if gene_id in bed_Dict.keys():
						gene = bed_Dict[gene_id]
					else:
						gene = line[2].strip()
					pvalue = line[3].strip()
					N = line[4].strip()
					rsid = line[5].strip()
					w.write(CHROM+"\t"+position+"\t"+gene+"\t"+pvalue+"\t"+N+"\t"+rsid+"\n")
			f.close()
			w.close()
if __name__ == '__main__':
	main(args)

