main(){
	global_parameters
	#create_folder
	#split_3aQTLs
	#run_sential_SNPs
	#run_coloc_analysis
	#sigSentinalSnp_all
	find_colocalizedGene

}


global_parameters(){
	CUR_DIR=`pwd`
	#origi_GWASs=$CUR_DIR/input/GWASs
	#origi_microArray=$CUR_DIR/input/microArray
	origi_transcript=$CUR_DIR/input/transcipt
	changed_GWASs=$CUR_DIR/output/GWASs
	#changed_microArray=$CUR_DIR/output/microArray
	changed_aQTLs=$CUR_DIR/output/aQTLs_transcript_hg19
	changed_transcript=$CUR_DIR/output/transcript
	tempdata=$CUR_DIR/output/tempdata
}


create_folder(){
	mkdir -p input/
	mkdir -p $tempdata
	mkdir -p $changed_GWASs
	mkdir -p $changed_microArray
	mkdir -p $changed_aQTLs
	mkdir -p $changed_transcript
	mkdir -p $CUR_DIR/output/analysis
	mkdir -p $CUR_DIR/results_transcript
}


split_3aQTLs(){
	mkdir -p $CUR_DIR/output/aQTLs_transcript
	mkdir -p $CUR_DIR/qsub
for file in $CUR_DIR/input/*.cis_aqtl_hg19.txt
	do
		tissueName=`echo "$file" |awk -F "/" '{print $NF;exit}'|awk -F".cis_aqtl_hg19" '{print $1;exit}'`
		cd $CUR_DIR/output/aQTLs_transcript
		if [[ ! -d "$tissueName" ]]
			then
				mkdir -p $CUR_DIR/output/aQTLs_transcript/${tissueName}
				echo '#!/bin/bash
#PBS -N split_'${tissueName}'
#PBS -q fat-1
#PBS -l nodes=1:ppn=16
#PBS -V
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
#module load R/3.6.2-anaconda3
cd $CUR_DIR/input/
head -n1 $file > $CUR_DIR/output/aQTLs_transcript/'$tissueName'_temptitle.txt
tail -n +2 $file |sort -k3 > $CUR_DIR/output/aQTLs_transcript/'$tissueName'_tempsorted.txt
wait
cat $CUR_DIR/output/aQTLs_transcript/'$tissueName'_temptitle.txt $CUR_DIR/output/aQTLs_transcript/'$tissueName'_tempsorted.txt >  $CUR_DIR/output/aQTLs_transcript/'$tissueName'_tempmerged.txt
wait
python2 $CUR_DIR/bin/extract_gene.py $CUR_DIR/output/aQTLs_transcript/'$tissueName'_tempmerged.txt $CUR_DIR/output/aQTLs_transcript/'$tissueName'
wait
rm  $CUR_DIR/output/aQTLs_transcript/'$tissueName'_temptitle.txt
rm  $CUR_DIR/output/aQTLs_transcript/'$tissueName'_tempsorted.txt
rm  $CUR_DIR/output/aQTLs_transcript/'$tissueName'_tempmerged.txt

wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
#module unload R/3.6.2-anaconda3
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait'>$CUR_DIR/qsub/split_${tissueName}.sh 
wait
cd $CUR_DIR/qsub/log
#qsub $CUR_DIR/qsub/split_${tissueName}.sh &
fi
echo $tissueName" is finished"
done
}



run_sential_SNPs(){
	for GWASfolder in $changed_GWASs/*
	do
		if [ "`ls $GWASfolder`" != "" ]
		then
			echo $GWASfolder
			echo $var
			trait=`echo "$GWASfolder" | awk -F"/" '{print $NF }'`
			echo ${trait}
			#mkdir -p $CUR_DIR/output/sentinalSNP/${trait}/src_pvalue_maf
			sentinalSNP_folder= $CUR_DIR/output/sentinelSNP/"$trait"
			cd $sentinalSNP_folder

			if [[ ! -f "sentinalSNP.txt" ]] && [[ ! -f "sentinalSNP_was_empty.txt" ]]
			then
				sentinalSNP_exist=0
				echo '#!/bin/bash
#PBS -N sen_'$trait'
#PBS -q cu-1
#PBS -l nodes=1:ppn=16
#PBS -V
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"


module load R/3.6.2-anaconda3
cp '$CUR_DIR'/src_pvalue_maf/get_sentinal_SNP.R '$CUR_DIR'/src_pvalue_maf/import.R '$sentinalSNP_folder'/src_pvalue_maf/
cd '$sentinalSNP_folder'/src_pvalue_maf/
Rscript get_sentinal_SNP.R --file '$GWASfolder'/GLGC_* &
wait
mv sentinalSNP*.txt '$sentinalSNP_folder'/ &
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

module unload intel_mpi/2018_u1
module unload R/3.6.2-anaconda3
wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
'> $sentinalSNP_folder/qsubSentinalSNP.sh &
				wait
				#qsub $sentinalSNP_folder/qsubSentinalSNP.sh &
				wait
			fi
		fi
	echo $trait" is all finished"
	#fi
	done


}

run_coloc_analysis(){
	var=1
	for GWASfolder in $changed_GWASs/*
	do
		if [ "`ls $GWASfolder`" != "" ]
		then
			echo $GWASfolder
			echo $var
			trait=`echo "$GWASfolder" | awk -F "/" '{print $NF }'`
			echo ${trait}
			#mkdir -p $CUR_DIR/output/sentinalSNP/${trait}/src_pvalue_maf
			sentinalSNP_folder= $CUR_DIR/output/sentinelSNP/"$trait"
			cd $sentinalSNP_folder

			if [[ ! -f "sentinalSNP.txt" ]] && [[ ! -f "sentinalSNP_was_empty.txt" ]]
			then
				sentinalSNP_exist=0
				echo '#!/bin/bash
#PBS -N sen_'$trait'
#PBS -q cu-1
#PBS -l nodes=1:ppn=16
#PBS -l walltime=100:00:00
#PBS -V
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

module load intel_mpi/2018_u1
module load R/3.6.2-anaconda3

cp '$CUR_DIR'/src_pvalue_maf/get_sentinal_SNP.R '$CUR_DIR'/src_pvalue_maf/import.R '$sentinalSNP_folder'/src_pvalue_maf/
cd '$sentinalSNP_folder'/src_pvalue_maf/
Rscript get_sentinal_SNP.R --file '$GWASfolder'/GLGC_*
wait
mv sentinalSNP*.txt '$sentinalSNP_folder'
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

module unload intel_mpi/2018_u1
module unload R/3.6.2-anaconda3
wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
'> $sentinalSNP_folder/qsubSentinalSNP.sh
				wait
				#qsub $sentinalSNP_folder/qsubSentinalSNP.sh
				wait
			else
				if [[ -f "sentinalSNP_was_empty.txt" ]]
				then
					sentinalSNP_exist=2
				elif [[ -f "sentinalSNP.txt" ]]
				then
					sentinalSNP_exist=1
				fi
			fi

			while (( sentinalSNP_exist == 0 ))
			do
				if [[ -f "sentinalSNP.txt" ]]
				then
					sentinalSNP_exist=1
				elif [[ -f "sentinalSNP_was_empty.txt" ]]
				then
					sentinalSNP_exist=2
				else
					sleep 1m
				fi
			done

			if [[ $sentinalSNP_exist == 2 ]]
			then
				echo $trait >> ../trait_No_sentinalSNP.txt
				continue
			fi
		fi


		if [ "`ls $GWASfolder`" != "" ]
		then
			trait=`echo "$GWASfolder" | awk -F"/" '{print $NF }'`
			#if [ $trait = "1210" ]
			#then

			for folder in $changed_aQTLs/*
			do
				if [ "`ls $folder`" != "" ]
				then
					tissue=`echo "$folder" | awk -F"/" '{print $NF }'`
			
					mkdir -p $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"
					if [ -d $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"/src_pvalue_maf/ ]
					then
						rm -r $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"/src_pvalue_maf/
						wait
					fi
					cp -r $CUR_DIR/src_pvalue_maf/ $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"
					wait
					cd $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"
					if [[ -f "summary_table.txt" ]] && [[ -f "coloc_bed_table.BED" ]]
					then
						echo -e "$tissue\t"$trait" exists"
						continue
					fi
					if [[ ! -f "summary_table.txt" ]]
					then
					cd $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"/src_pvalue_maf/
					echo '#!/bin/bash
#PBS -N '$trait'_'$tissue'
#PBS -q cu-1
#PBS -l nodes=1:ppn=16
#PBS -l walltime=100:00:00
#PBS -V
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

module load intel_mpi/2018_u1
module load R/3.6.2-anaconda3

cd '$CUR_DIR'/output/results_transcript/'$trait'/'$tissue'_'$trait'/src_pvalue_maf/
Rscript wrapper.R ./ '$CUR_DIR'/output/aQTLs_transcript_hg19/'$tissue'/ '$CUR_DIR'/output/GWASs/'$trait'/ '$sentinalSNP_folder'/sentinalSNP.txt

wait
mv summary_table.txt coloc_bed_table.BED ../
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

module unload intel_mpi/2018_u1
module unload R/3.6.2-anaconda3
wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
'> $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"/qsubWrapper.sh
					wait
					cd $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"
					qsub $CUR_DIR/output/results_transcript/"$trait"/"$tissue"_"$trait"/qsubWrapper.sh &
					wait
					fi
				fi
		process=`qstat -u |wc -l`
		while (( process >= 300))
		do
			echo "Current number of jobs is larger than 300"
			echo "Wait another 10 minutes"
			sleep 10m
			process=`qstat -u |wc -l`
		done
		if [ $var -lt 1 ]
		then
			let "var+=1"
			continue
		fi
		if [ $var -gt 100000000 ]
		then
			break
		fi
		let "var += 1"			

			done
			#fi
			echo $trait" is all finished"
		fi

	done



}

sigSentinalSnp_all(){
echo '#!/bin/bash
#PBS -N sigSentinalSnp
#PBS -q cu-1
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -l walltime=100:00:00
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"        


cd '$CUR_DIR'/bin
        #Use this script to calculate the ratio of significant sentinal snp in all tissue and each trait
        python sigSentinalSnp_all_105.py
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/' > $CUR_DIR/qsub/sigSentinalSnp_1e5.sh
			 cd $CUR_DIR/qsub
                         #qsub sigSentinalSnp.sh
}




find_colocalizedGene(){
echo '#!/bin/bash
#SBATCH --job-name='find_colocalizedGene'
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=%j.err
#SBATCH --output=%j.out

##################################

I=`whoami`
CURDIR=`pwd`
NODES=`scontrol show hostname $SLURM_JOB_NODELIST`

for i in $NODES

do

echo "$i:$SLURM_NTASKS_PER_NODE" >> $CURDIR/nodelist.$SLURM_JOB_ID

done

echo $SLURM_NPRCOS

echo "process will start at : "

date

echo "++++++++++++++++++++++++++++++++++++++++"
	cd '$CUR_DIR'/bin
	#python find_colocalizedGene_v2.py ## extract the sig terms
	python find_colocalizedGene_v3_all.py ## extract all terms
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/' > $CUR_DIR/qsub/find_colocalizedGene.sh
                        cd $CUR_DIR/qsub 
			sbatch find_colocalizedGene.sh
}


main
