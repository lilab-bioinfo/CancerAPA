#!/bin/python
#Use this script to calculate the proportion of significant sentinal snp.
#The default run path is You-Path-To-Project/bin/


from __future__ import division
import sys
import os

file_results=os.getcwd()+"/../output/results_transcript"
trait_file=os.listdir(file_results)
file_trait=sorted(trait_file)
file_sentinalSNP=os.getcwd()+"/../output/sentinalSNP"
path_output=file_results+"/../analysis/allTissueProportion_SigSentinalSnp_range1e6_hyp4_0.75.txt"
file_output=open(path_output,"w+")

#Find all the trait type and for each trait, set a list variable
traits=[]
for each_trait in file_trait:
	path_tissue=file_results+"/"+each_trait
	trait_folder=os.listdir(path_tissue)
	for each_trait_tissue in trait_folder:
		tissueName='_'.join(each_trait_tissue.split('_')[:-1])
		file_sentinal=file_sentinalSNP+"/"+each_trait+"/sentinalSNP.txt"
		if not os.path.exists(file_sentinal):
			continue
		trait=each_trait_tissue.split('_')[-1]
		if not trait in traits:
			traits.append(trait)
			locals()[trait]=[]
			file2=open(file_sentinal,"r")
			f2=file2.readlines()
			file2.close()
			rangeSig=1e5
			for line in f2:
				line_split=line.split('\t')
				if line_split[1].isdigit():
					sentinalSnpTSS=line_split[1]
					snpChr=line_split[0][1:-1]
					bottom=int(sentinalSnpTSS)-rangeSig
					top=int(sentinalSnpTSS)+rangeSig
					flag=0
					snp=line_split[3][1:-1]
					sentinalSnp=[bottom,top,flag,snpChr,snp]
					print(sentinalSnp)
					locals()[trait].append(sentinalSnp)

for each_trait in file_trait:
	print(each_trait)
	path_tissue=file_results+"/"+each_trait
	trait_folder=os.listdir(path_tissue)
	for each_trait_tissue in trait_folder:
		tissueName='_'.join(each_trait_tissue.split('_')[:-1])
		print(tissueName)
		trait=each_trait_tissue.split('_')[-1]
		file_summary=path_tissue+"/"+each_trait_tissue+"/summary_table.txt"
		file_coloc=path_tissue+"/"+each_trait_tissue+"/coloc_bed_table.BED"
		print(file_summary)
		print(file_coloc)
		if not os.path.exists(file_summary):
			continue
		if not os.path.exists(file_coloc):
			continue
		if trait not in traits:
			continue
		file1=open(file_summary,"r")
		file2=open(file_coloc,"r")
		f1=file1.readlines()
		f2=file2.readlines()
		file1.close()
		file2.close()
		geneChrDic={}
		for line in f2:
			line_split=line.split('\t')
			if len(line_split)>3:
				key=line_split[3]
				value=line_split[0]
				if key not in geneChrDic:
					geneChrDic[key] = value
		GeneTSS=[]
		for line in f1:
			line_split=line.split('\t')
			if line_split[5] !="NA" and line_split[6] !="NA":
				hyp3=float(line_split[5])
				hyp4=float(line_split[6])
				#if hyp4>0.75:
				#if hyp4>0.75:#by HYM
				if hyp4>0.75 and hyp4/(hyp3+hyp4)>=0.9: 
					geneName=line_split[0].split('_')[0]
					if geneName != "summary" and geneName != "me" and geneName in geneChrDic:
						geneChr=geneChrDic[geneName]
						geneTSS=line_split[0].split('_')[2]
						GeneTSS.append([geneChr,geneTSS])
		for i in range(len(locals()[trait])):
			if locals()[trait][i][2] ==1:
				continue
			else:
				for gene in GeneTSS:
					try:
						if gene[0] == locals()[trait][i][3] and int(gene[1]) > locals()[trait][i][0] and int(gene[1]) < locals()[trait][i][1]:
							locals()[trait][i][2] =1
							GeneTSS.remove(gene)
							break
					except Exception:
						continue

for trait in traits:
	totalSentinalSnp=len(locals()[trait])
	sigSentinalSnpNum=0
	colocsnp=[]
	for eachTrait in locals()[trait]:
		if eachTrait[2] == 1:
			sigSentinalSnpNum += 1
			colocsnp.append(eachTrait[4])
	#print locals()[trait]
	#proportionSigSentinalSnp = sigSentinalSnpNum/totalSentinalSnp
	file_output.write(trait+"\t"+str(sigSentinalSnpNum)+"\t"+str(totalSentinalSnp)+"\t"+','.join(colocsnp)+"\n") 


file_output.close()
