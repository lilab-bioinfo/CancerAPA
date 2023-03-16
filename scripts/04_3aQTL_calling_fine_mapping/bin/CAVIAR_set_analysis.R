library(data.table)
library(dplyr)

options(echo=TRUE) # if you want to see commands in output file
#args <- commandArgs(TRUE)

#tissue<-args[1]
#print(tissue)
#Transcript<-args[2]
#print(Transcript)

list1<-fread("/lustre/home/hchen/2021-10-31-cancer-GWAS/2022-09-19-CAVIAR-aQTL-finemapping/input/GTEx_gene_list_109117_149117.txt",sep="\t",header=F)
list2<-fread("/lustre/home/hchen/2021-10-31-cancer-GWAS/2022-09-19-CAVIAR-aQTL-finemapping/input/GTEx_gene_list_1_9999.txt",sep="\t",header=F)
list3<-fread("/lustre/home/hchen/2021-10-31-cancer-GWAS/2022-09-19-CAVIAR-aQTL-finemapping/input/GTEx_gene_list_96543_100000.txt",sep="\t",header=F)
genelist<-rbind(list1,list2,list3)


#set=c()
for (i in 1:dim(genelist)[1]){
	tissue=genelist[i,1]
	print(tissue)
	Transcript=genelist[i,2]
	print(Transcript)
	outputdir="/lustre/home/hchen/2021-10-31-cancer-GWAS/2022-09-19-CAVIAR-aQTL-finemapping/output/"
	outputfile=paste0("/lustre/home/hchen/2021-10-31-cancer-GWAS/2022-09-19-CAVIAR-aQTL-finemapping/output_sets/hchen/",tissue,"_",Transcript,"_3aQTL.cts_95post")
	if(file.exists(outputfile)){
		next()
	}


	cts_post_file=paste0(outputdir,tissue,"/",Transcript,"/caviar/",tissue,"_3aQTL.cts_post") ### 99% credible sets with all SNPs
	cts_set_file=paste0(outputdir,tissue,"/",Transcript,"/caviar/",tissue,"_3aQTL.cts_set")  ## 95% credible sets

	if(!file.exists(cts_post_file) || !file.exists(cts_set_file) || file.info(cts_set_file)$size<2){
                next()
        }
	
	post_dat<-fread(cts_post_file,header=T,sep="\t")
	set_dat<-fread(cts_set_file,header=F)
	post_dat%>%filter(post_dat$SNP_ID %in% set_dat$V1) -> set_post_dat ## filter the posterior of 95% credible sets

	set_post_dat$tissue<-tissue ## add corresponding tissue and gene information
	set_post_dat$gene_name<-Transcript
	#outputfile=paste0("/lustre/home/hchen/2021-10-31-cancer-GWAS/2022-09-19-CAVIAR-aQTL-finemapping/output_sets/hchen/",tissue,"_",Transcript,"_3aQTL.cts_95post")
	#if(!file.exists(outputfile)){
	fwrite(set_post_dat,outputfile,quote=F,sep="\t")
	#}else{
	#	next()
	#}	
	#set<-rbind(set,set_post_dat)
}

#outputfile=paste0("/lustre/home/hchen/2021-10-31-cancer-GWAS/2022-09-19-CAVIAR-aQTL-finemapping/output_sets/GTEx_gene_1_9999_3aQTL.cts_95post")
#fwrite(set_post_dat,outputfile,quote=F,sep="\t")

