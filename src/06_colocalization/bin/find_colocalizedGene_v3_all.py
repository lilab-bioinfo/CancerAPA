#!/bin/python

#Find the colocalized gene in summary_table.txt of each tissue-trait with the standard of hypp4>0.75 and hypp4/(hypp3+hypp4) >0.9
#The default run path is You-Path-To-Project/bin/

import sys
import os
#from pathlib import Path

file_results=os.getcwd()+"/../output/results_transcript"
file_GWASs=os.getcwd()+"/../output/results_transcript"

file_results_tissue=os.listdir(file_results)
file_results_trait=os.listdir(file_GWASs)

file_results_trait=sorted(file_results_trait)

path_output=file_results+"/../analysis/All_colocalizedGene_3aQTLs.txt"

file_output=open(path_output,"w+")
file_output.write("GeneName"+"\t"+"SNP_Num"+"\t"+"HyPP0"+"\t"+"Hypp1"+"\t"+"Hypp2"+"\t"+"Hypp3"+"\t"+"Hypp4"+"\t"+"Tissue"+"\t"+"Trait"+"\n")


traitStart=1
traitFinished=80
#Analyze traits 
for i in range(traitStart-1,traitFinished):
	trait=file_results_trait[i]
	print(trait)
	my_file = Path(file_results+"/"+trait)
	if my_file.exists():
		file_results_tissue=os.listdir(file_results+"/"+trait)
		for tissueTrait in file_results_tissue:
			print(tissueTrait)
			tissueName=tissueTrait.split('_')[0]
			if len(tissueTrait.split('_')) > 2:
				for i in range(1,len(tissueTrait.split('_'))-1):
					tissueName=tissueName+"_"+tissueTrait.split('_')[i]
			traitName=tissueTrait.split('_')[-1]
			print(traitName)
			if traitName != trait:
				continue
			if not os.path.exists(file_results+"/"+trait+"/"+tissueTrait+"/summary_table.txt"):
				continue
			fileSummary=open(file_results+"/"+trait+"/"+tissueTrait+"/summary_table.txt","r")
			f1=fileSummary.readlines()
			fileSummary.close()
			for line in f1:
				line_split=line.split('\t')
				if line_split[5]!= "NA" and line_split[6]!="NA":
					hyp3=float(line_split[5])
					hyp4=float(line_split[6][:-1])
					if hyp4 >0:
						geneName=line_split[0].split('_')[0]
						key=tissueName+"_"+traitName+"_"+geneName
	#			if key in gene_Dic:
	#				value=gene_Dic[key]
	#				file_output.write(geneName+"\t"+line_split[1]+"\t"+line_split[2]+"\t"+line_split[3]+"\t"+line_split[4]+"\t"+line_split[5]+"\t"+line_split[6][:-1]+"\t"+tissueName+"\t"+traitName+"\t"+"Colocalized"+"\t"+value+"\n")
	#			else:
	#				value="NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"
	#				file_output.write(geneName+"\t"+line_split[1]+"\t"+line_split[2]+"\t"+line_split[3]+"\t"+line_split[4]+"\t"+line_split[5]+"\t"+line_split[6][:-1]+"\t"+tissueName+"\t"+traitName+"\t"+"No-Colocalized"+"\t"+value+"\n")
						file_output.write(geneName+"\t"+line_split[1]+"\t"+line_split[2]+"\t"+line_split[3]+"\t"+line_split[4]+"\t"+line_split[5]+"\t"+line_split[6][:-1]+"\t"+tissueName+"\t"+traitName+"\n")


file_output.close()
