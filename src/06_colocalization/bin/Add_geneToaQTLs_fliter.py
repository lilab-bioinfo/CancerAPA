# -*- coding: utf-8 -*-
#python ~/aQTL_coloc/bin/Add_geneToaQTLs_fliter.py -b ~/aQTL_coloc/bin/hg38_Refseq_id_from_UCSC.txt -p ~/aQTL_coloc/input/20210718_aQTLs_postionTohg19_for_wenyan -o ~/aQTL_coloc/input/20210728_aQTLs_postionTohg19_for_coloc 

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
	parser.add_argument("-b","--bedfile",help="hg38_Refseq_id_from_UCSC.txt pathfile")
	parser.add_argument("-p","--pathaqtl",help="aQTL pathfile")
	parser.add_argument("-o","--output",help="output pathfile")
	args = parser.parse_args()
	path_aqtl = args.pathaqtl
	bedfile = args.bedfile
	outputfile = args.output
	bed_Dict = {}
	files_aqtl = os.listdir(path_aqtl)
	bed = open(bedfile)
	for line in bed:
		line = line.split("\t")
		if line[0].strip()!="#name":
			gene_id = line[0].strip()
			bed_Dict[gene_id] = line[1].split(".")[0].strip()
	for file in files_aqtl:
		if not os.path.isdir(file):
			f = open(path_aqtl+"/"+file)
			w = open(outputfile+"/"+file,"w")
			w.write("chr"+"\t"+"position"+"\t"+"gene"+"\t"+"p-value"+"\t"+"N"+"\t"+"rsid"+"\n")
			for line in f:
				line = line.split("\t")
				if line[0].strip()!="gene":
					CHROM = line[1].strip()
					position = line[2].strip()
					gene_id = line[0].split("|")[0].split(".")[0].strip()
					if gene_id in bed_Dict.keys():
						gene = bed_Dict[gene_id]
					else:
						gene = line[0].split("|")[1].strip()
					pvalue = line[4].strip()
					N = line[7].strip()
					rsid = line[3].strip()
					w.write(CHROM+"\t"+position+"\t"+gene+"\t"+pvalue+"\t"+N+"\t"+rsid+"\n")
			f.close()
			w.close()
if __name__ == '__main__':
	main(args)

