#!/bin/python3

#python3 /media/disk2/eQTL/bin/gtf_get_bed1.py /media/disk2/eQTL/input/transcript/gencode.v31.annotation.gtf /media/disk2/eQTL/output/transcript/gencode.v31.annotation_gene.bed &

#Use this script to get the .bed file which contains:
	#   V1   | V2  | V3  | V4               | V5        | V6
	#   chr# | TSS | TES | common_gene_name | ensemblID | strand

import sys
from gtfparse import read_gtf
import numpy as np
import pandas as pd

path_input=sys.argv[1]
path_output=sys.argv[2]
file_input=read_gtf(path_input)
# filter DataFrame to gene
df_genes = file_input[file_input["feature"] == "gene"]
#file_input.close()
file_output=open(path_output,"w+")
output = df_genes[["seqname","start","end","gene_name","gene_id","strand"]]
output.to_csv(file_output, header=False,sep='\t',index=0 )	
#file_output.close()