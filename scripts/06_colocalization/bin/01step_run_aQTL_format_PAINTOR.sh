	var=1
	mkdir -p /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL/qsub/aQTL_format_PAINTOR
	mkdir -p /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL/output/aQTL_format_PAINTOR
	cd /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL
	#rm *_split*
	for file in /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL/input/output_aQTL/*.cis_eqtl_all_approachb.txt.gz
	do
		fileName=`echo "$file" |awk -F"/" '{print $NF;exit}'|awk -F"." '{print $1;exit}'`
		echo '#!/bin/bash
#PBS -N '$fileName'
#PBS -q fat-1
#PBS -l nodes=1:ppn=2
#PBS -l walltime=200:00:00
#PBS -V
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$
module load R/3.6.2-anaconda3
module load samtools/1.10-gcc-4.8.5
module load htslib/1.10.5-gcc-4.8.5
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
module load R/3.6.2-anaconda3
module load openmpi/openmpi-3.1.2-intel-2019-update5
module load gcc/9.2.0
cd /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL
python /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL/bin/aQTL_format_PAINTOR_hg38.py -c '$file' -i /lustre/home/ymhu/TWAS/input/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz -o /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL

wait
echo "processs will sleep 5s"
sleep 5
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

'> /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL/qsub/aQTL_format_PAINTOR/${fileName}.sh
		#wait
	qsub /lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL/qsub/aQTL_format_PAINTOR/${fileName}.sh
	wait
					#fi
		echo $fileName" is all finished"
		process=`qstat | grep "hchen"|wc -l`
		while (( process >= 590))
		do
			echo "Current number of jobs is larger than 590"
			echo "Wait another 5 minutes"
			sleep 5m
			process=`qstat | grep "hchen"|wc -l`
		done
		if [ $var -lt 1 ]
		then
			let "var+=1"
			continue
		fi
		if [ $var -gt 1000000000 ]
		then
			break
		fi
		let "var += 1"

		done
