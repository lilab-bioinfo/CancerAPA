#!/bin/bash

# --------------------- Functions ----------------------
# main function

main(){
#munge_sumstats
#ldsc_annot
#paste_all
partitioned_heritability_allgwas_baseline
}

munge_sumstats(){
usr=`pwd`
usr1=`echo $usr|awk -F "/" '{print $4}'`
gwas=/lustre/home/${usr1}/GWAS_SMR_format
path=/lustre/home/${usr1}/Partitioned-Heritability	
mkdir -p ${path}/output
mkdir -p ${path}/script/munge_sumstats
mkdir -p ${path}/input/munge_sumstats
ldsc=/lustre/home/${usr1}/src/ldsc-master

for file in ${gwas}/*.txt;do
		tissueName=`echo "$file" |awk -F"/" '{print $NF;exit}'|awk -F".txt" '{print $1;exit}'`
		cd ${path}/input/munge_sumstats
		if [[ ! -f ${tissueName}.sumstats.gz ]]
			then	
	echo '#!/bin/bash
#SBATCH --job-name='${tissueName}'
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --error='${tissueName}'.err
#SBATCH --output='${tissueName}'.out
##########################################

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

python '${ldsc}'/munge_sumstats.py --sumstats '${gwas}'/'${tissueName}'.txt --merge-alleles '${ldsc}'/test/w_hm3.snplist --out '${path}'/input/munge_sumstats/'${tissueName}'
 
wait
echo "processs will sleep 30s"
sleep 3
echo "process end at : "
date

rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait'> ${path}/script/munge_sumstats/${tissueName}.sh 
wait
cd ${path}/script/munge_sumstats
sbatch ${path}/script/munge_sumstats/${tissueName}.sh &
echo $tissueName" is finished"
  process=`squeue | grep "${usr1}"|wc -l`
		while (( process >= 99))
		do
			echo "Current number of jobs is larger than 100"
			echo "Wait another 1 minutes"
			sleep 10
			process=`squeue | grep "${usr1}"|wc -l`
		done
fi		
done

}

ldsc_annot(){
usr=`pwd`
usr1=`echo $usr|awk -F "/" '{print $4}'`
path=/lustre/home/${usr1}/Partitioned-Heritability
mkdir -p ${path}/script/ldsc_annot
mkdir -p ${path}/input/ldsc_annot/aQTL
base_dirt_aqtl=${path}/input/add_rsid/
geno_dirt=/lustre/home/${usr1}/data/ref/anno/hg38/plink_files/1000G.EUR.hg38.
result=${path}/input/ldsc_annot/aQTL
script=${path}/script/ldsc_annot

for file in ${base_dirt_aqtl}/*.rs;do
		tissueName=`echo "$file" |awk -F"/" '{print $NF;exit}'|awk -F".rs" '{print $1;exit}'`
for chr in $(seq 1 22); do		
echo '
library(data.table)
aQTL=fread(paste0("'${base_dirt_aqtl}'","'${tissueName}'.rs"),head=F,stringsAsFactors=F,data.table=F)
aQTL1=na.omit(unique(aQTL[,1]))
bim=fread(paste0("'${geno_dirt}'","'${chr}'",".bim"),head=F,stringsAsFactors=F,data.table=F)
bim$aQTL1=0; 
index3=match(aQTL1,bim$V2,nomatch=0)
bim$aQTL1[index3]=1;
anot_aQTL=bim[,c(1,4,2,3,7)]
colnames(anot_aQTL)=c("CHR","BP","SNP","CM","'${tissueName}'")
write.table(anot_aQTL[,5],paste0("'${result}'/'${tissueName}'.","'${chr}'",".annot"),row=F,col=T,quo=F,sep="\t")

'> ${script}/${tissueName}.${chr}.r
		
				echo '#!/bin/bash
#SBATCH --job-name='${tissueName}'.'${chr}'
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --error='${tissueName}'.'${chr}'.err
#SBATCH --output='${tissueName}'.'${chr}'.out
##########################################

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

module load R/3.6.2-anaconda3

cd '${script}'
R CMD BATCH --no-save '${script}'/'${tissueName}'.'${chr}'.r

wait
sed -i '\''1c '${tissueName}''\'' '${result}'/'${tissueName}'.'${chr}'.annot

wait
echo "processs will sleep 30s"
sleep 3
echo "process end at : "
date
module unload R/3.6.2-anaconda3
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait'> ${script}/${tissueName}.${chr}.sh 
wait
cd ${script}
sbatch ${script}/${tissueName}.${chr}.sh &
echo ${tissueName}.${chr}" is finished"
  process=`squeue | grep "${usr1}"|wc -l`
		while (( process >= 99))
		do
			echo "Current number of jobs is larger than 100"
			echo "Wait another 1 minutes"
			sleep 10
			process=`squeue | grep "${usr1}"|wc -l`
		done
done
done
}


paste_all(){
usr=`pwd`
usr1=`echo $usr|awk -F "/" '{print $4}'`
path=/lustre/home/${usr1}/Partitioned-Heritability
mkdir -p ${path}/input/paste_all/aQTL
mkdir -p ${path}/script/paste_all
result=${path}/input/paste_all/aQTL
script=${path}/script/paste_all
bin=/lustre/home/${usr1}/src/ldsc-master
geno_dirt=/lustre/home/${usr1}/data/ref/anno/hg38/plink_files/1000G.EUR.hg38

for i in $(seq 1 22); do
echo '#!/bin/bash
#SBATCH --job-name='${i}'
#SBATCH --partition=fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --error='${i}'.err
#SBATCH --output='${i}'.out
##########################################

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

#for 3aQTL
#paste '${path}'/input/ldsc_annot/aQTL/*.'${i}'.annot >> '${result}'/aQTL.'${i}'.annot
#gzip '${result}'/aQTL.'${i}'.annot

source activate /lustre/home/${usr1}/anaconda3/envs/ldsc
python '${bin}'/ldsc.py --l2 --bfile '${geno_dirt}'.'${i}' --ld-wind-cm 1 --annot '${result}'/aQTL.'${i}'.annot.gz --thin-annot --out '${result}'/aQTL.'${i}' --print-snps '${bin}'/list.txt

conda deactivate 
wait
echo "processs will sleep 30s"
sleep 3
echo "process end at : "
date

wait' > ${script}/${i}.sh

  wait
  cd ${script}
  sbatch ${script}/${i}.sh
  echo ${i}" is all finished"
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

partitioned_heritability_allgwas_baseline(){
usr=`pwd`
usr1=`echo $usr|awk -F "/" '{print $4}'`
path=/lustre/home/${usr1}/Partitioned-Heritability
mkdir -p ${path}/script/partitioned_heritability_allgwas_baseline
mkdir -p ${path}/output
script=${path}/script/partitioned_heritability_allgwas_baseline
result=${path}/output
aqtl=${path}/input/paste_all/aQTL
update=/lustre/home/${usr1}/Partitioned-Heritability/input/munge_sumstats
envs=/lustre/home/${usr1}/anaconda3/envs
bin=/lustre/home/${usr1}/src/ldsc-master
ref=/lustre/home/${usr1}/data/ref/anno/hg38

for file in ${update}/*.sumstats.gz;do
		tissueName=`echo "$file" |awk -F"/" '{print $NF;exit}'|awk -F".sumstats.gz" '{print $1;exit}'`
		cd ${path}/output
		if [[ ! -f ${line1}.results ]]
		then
echo '#!/bin/bash
#SBATCH --job-name='${tissueName}'
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --error='${tissueName}'.err
#SBATCH --output='${tissueName}'.out
##########################################

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

source activate '${envs}'/ldsc

python '${bin}'/ldsc.py --h2 '${update}'/'${tissueName}'.sumstats.gz --ref-ld-chr '${ref}'/baseline_v1.2/baseline.,'${aqtl}'/aQTL. --w-ld-chr '${ref}'/weights/weights.hm3_noMHC. --overlap-annot --print-coefficients --frqfile-chr '${ref}'/plink_files/1000G.EUR.hg38. --out '${result}'/'${tissueName}'

conda deactivate 

wait
echo "processs will sleep 30s"
sleep 5
echo "process end at : "
date

wait' > ${script}/${tissueName}.sh

  wait
  cd ${script}
  sbatch ${script}/${tissueName}.sh
  echo ${tissueName}" is all finished"
  fi
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




main
