#!/bin/bash

# --------------------- Functions ----------------------
# main function


main(){
#generate_anno
fgwas
}

generate_anno(){
usr=`pwd`
usr1=`echo $usr|awk -F "/" '{print $4}'`
path=/lustre/home/${usr1}/fgwas-cancer
mkdir -p ${path}/output
mkdir -p ${path}/script/generate_anno
mkdir -p ${path}/input/generate_anno
gwas=/lustre/home/${usr1}/fgwas-cancer/coloc_hg38
eqtl=/lustre/home/${usr1}/fgwas-cancer/eqtl_sig
sqtl=/lustre/home/${usr1}/fgwas-cancer/sqtl_sig
aqtl=/lustre/home/${usr1}/fgwas-cancer/aqtl_sig

for GWASfolder in ${gwas}/CG*;do
		tissueName=`echo "$GWASfolder" | awk -F"/" '{print $NF }'`
echo '
library("dplyr")
library(stringr)
eqtl <- read.table("'${gwas}'/'${tissueName}'/GLGC_'${tissueName}'_result.txt", header = T,sep="\t",check.names=F)
eqtl$Z=eqtl$beta/eqtl$se
eqtl1=eqtl[,c("rsid","chr","position","maf","Z","N")]
colnames(eqtl1)=c("SNPID","CHR","POS","F","Z","N")
eqtl4=eqtl1[order(as.numeric(sub("\\chr+", "", eqtl1$CHR)),eqtl1$POS),]

filePath <- list.files(path = "'${eqtl}'",pattern = "*.snp.uniq", full.names=T)
for(i in 1:length(filePath)){
g_file <- filePath[i]
data <- read.table(g_file,header=F) %>% as.data.frame %>% mutate(V2="1")
fileNames <- basename(g_file)
colnames(data)=c("SNPID",paste0(as.character(strsplit(as.character(fileNames), split = "[.]")[[1]][1]),"_eQTL"))
eqtl4 <- eqtl4 %>% left_join(data,by=c("SNPID"="SNPID"))
eqtl4[is.na(eqtl4)] <- 0
}

filePath <- list.files(path = "'${sqtl}'",pattern = "*.snp.uniq", full.names=T)
for(i in 1:length(filePath)){
g_file <- filePath[i]
data <- read.table(g_file,header=F) %>% as.data.frame %>% mutate(V2="1")
fileNames <- basename(g_file)
colnames(data)=c("SNPID",paste0(as.character(strsplit(as.character(fileNames), split = "[.]")[[1]][1]),"_sQTL"))
eqtl4 <- eqtl4 %>% left_join(data,by=c("SNPID"="SNPID"))
eqtl4[is.na(eqtl4)] <- 0
}

filePath <- list.files(path = "'${aqtl}'",pattern = "*.snp.uniq", full.names=T)
for(i in 1:length(filePath)){
g_file <- filePath[i]
data <- read.table(g_file,header=F) %>% as.data.frame %>% mutate(V2="1")
fileNames <- basename(g_file)
colnames(data)=c("SNPID",paste0(as.character(strsplit(as.character(fileNames), split = "[.]")[[1]][1]),"_aQTL"))
eqtl4 <- eqtl4 %>% left_join(data,by=c("SNPID"="SNPID"))
eqtl4[is.na(eqtl4)] <- 0
}

write.table(eqtl4, "'${path}'/input/generate_anno/'${tissueName}'.txt",sep =" ", row.names =FALSE, col.names =TRUE, quote =FALSE)


'> ${path}/script/generate_anno/${tissueName}.r
		
				echo '#!/bin/bash
#SBATCH --job-name='${tissueName}'
#SBATCH --partition=fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error='${tissueName}'.err
#SBATCH --output='${tissueName}'.out
##########################################

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

cd '${path}'/script/generate_anno
R CMD BATCH --no-save '${path}'/script/generate_anno/'${tissueName}'.r
gzip '${path}'/input/generate_anno/'${tissueName}'.txt
 
wait
echo "processs will sleep 30s"
sleep 3
echo "process end at : "
date
#module unload R/3.6.2-anaconda3
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait'> ${path}/script/generate_anno/${tissueName}.sh 
wait
cd ${path}/script/generate_anno
sbatch ${path}/script/generate_anno/${tissueName}.sh &
echo $tissueName" is finished"
  process=`squeue | grep "${usr1}"|wc -l`
		while (( process >= 99))
		do
			echo "Current number of jobs is larger than 100"
			echo "Wait another 1 minutes"
			sleep 10
			process=`squeue | grep "${usr1}"|wc -l`
		done
done

}


fgwas(){
usr=`pwd`
usr1=`echo $usr|awk -F "/" '{print $4}'`
path=/lustre/home/${usr1}/fgwas-cancer
mkdir -p ${path}/output
mkdir -p ${path}/script/fgwas
for file in ${path}/input/generate_anno/*.txt.gz;do
		tissueName=`echo "$file" |awk -F"/" '{print $NF;exit}'|awk -F".txt.gz" '{print $1;exit}'`
while read line
do 		
cd ${path}/output
		if [[ ! -f ${line}_${tissueName}.params ]]
			then	
				echo '#!/bin/bash
#SBATCH --job-name='${line}'_'${tissueName}'
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --error='${line}'_'${tissueName}'.err
#SBATCH --output='${line}'_'${tissueName}'.out
##########################################

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

fgwas -i '${path}'/input/generate_anno/'${tissueName}'.txt.gz -w '${line}' -o '${path}'/output/'${line}'_'${tissueName}'
 
wait
echo "processs will sleep 30s"
sleep 3
echo "process end at : "
date
#module unload R/3.6.2-anaconda3
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait'> ${path}/script/fgwas/${line}_${tissueName}.sh 
wait
cd ${path}/script/fgwas
sbatch ${path}/script/fgwas/${line}_${tissueName}.sh &
echo ${line}_$tissueName" is finished"
  process=`squeue | grep "${usr1}"|wc -l`
		while (( process >= 99))
		do
			echo "Current number of jobs is larger than 100"
			echo "Wait another 1 minutes"
			sleep 10
			process=`squeue | grep "${usr1}"|wc -l`
		done
fi		
done < anno.list
done

}


main
