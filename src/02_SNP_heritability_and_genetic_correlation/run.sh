main(){
      global_parameters
      creat_folder
      run_menge_sumstats
      run_ldsc_heritability
      run_ldsc_genetic_correlation
}

global_parameters(){
	CUR_DIR=$(cd `dirname $0` && pwd -P)
	LDSC_GWASs=$CUR_DIR/input
	LDSC_output=$CUR_DIR/result
}


create_folder(){
	mkdir -p input/
	mkdir -p $LDSC_GWASs
	mkdir -p $LDSC_output
	mkdir -p $CUR_DIR/software/ldsc/qusb
	mkdir -p $LDSC_output/ldsc_pair
}

run_munge_sumstats(){
	var=1
	mkdir -p $LDSC_output/20210729_munge_sumstats
	for file in {CG0051,CG0052,CG0102,CG0251,CG0255,CG0256,CG0258,CG0276,CG0361,CG0006,CG0007,CG0054,CG0083,CG0103,CG0104,CG0123,CG0124,CG0118,CG0133,CG0134,CG0154,CG0155,CG0259,CG0265,CG0266,CG0277,CG0055,CG0056,CG0116,CG0117,CG0057,CG0003,CG0312,CG0008,CG0012,CG0027,CG0030,CG0033,CG0042,CG0058,CG0068,CG0071,CG0084,CG0085,CG0089,CG0105,CG0125,CG0128,CG0119,CG0356,CG0308,CG0139,CG0142,CG0150,CG0156,CG0171,CG0172,CG0175,CG0178,CG0181,CG0184,CG0198,CG0201,CG0204,CG0359,CG0218,CG0237,CG0366,CG0260,CG0281,CG0291,CG0294,CG0301,CG0320,CG0323,CG0353,CG0038,CG0317,CG0362,CG0419,CG0011,CG0398,CG0388,CG0418,CG0077,CG0079,CG0080,CG0108-B,CG0191-B,CG0193-B,CG0381,CG0086,CG0087,CG0417,CG0435,CG0109,CG0127,CG0385,CG0130,CG0416,CG0421,CG0431,CG0432,CG0370,CG0371,CG0391,CG0357,CG0423,CG0136,CG0137,CG0386,CG0212,CG0239,CG0425,CG0436,CG0437,CG0146,CG0192,CG0194,CG0196,CG0197,CG0394,CG0397,CG0209,CG0422,CG0207,CG0208,CG0210,CG0211,CG0213,CG0214,CG0215,CG0424,CG0216,CG0360,CG0221,CG0222,CG0223,CG0224,CG0225,CG0226,CG0227,CG0228,CG0392,CG0393,CG0428,CG0429,CG0238,CG0387,CG0261,CG0268,CG0300,CG0304,CG0420,CG0327,CG0328,CG0342,CG0343,CG0344,CG0382,CG0383,CG0384,CG0426,CG0427,CG0434,CG0039,CG0345,CG0346,CG0364}
	do
		GWASfile=$LDSC_GWASs/GLGC_${file}_result.txt
	
		if [[ -f "$GWASfile" ]]
		then
			echo $GWASfile
			
			trait=`echo "$GWASfile" | awk -F "/" '{print $NF }'`
			cd $LDSC_output/20210729_munge_sumstats
			CGID=`echo "$trait" | awk -F "_" '{print $2 }'`

			if [[ ! -f "$trait".sumstats.gz ]]
			then
				echo '
#!/bin/bash
#PBS -N munge_'$CGID'
#PBS -q cu-1
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
#source activate ldsc
module load R/3.6.2-anaconda3
source '$CUR_DIR'/.bashrc
source activate '$CUR_DIR'/.conda/envs/ldsc
cd '$LDSC_output'/20210729_munge_sumstats
'$CUR_DIR'/software/ldsc/munge_sumstats.py --sumstats '$GWASfile' --out '$trait' --a1-inc --merge-alleles '$CUR_DIR'/software/ldsc/w_hm3.snplist --chunksize 10000000000

wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

'> $CUR_DIR/software/ldsc/qusb/"$trait".sh
					wait
					qsub $CUR_DIR/software/ldsc/qusb/"$trait".sh
					wait
					#fi
				fi
			
			echo $trait" is all finished"
		fi

done

}


run_ldsc_heritability(){
	mkdir -p $LDSC_output/ldsc_heritability
	for GWASfile in $LDSC_output/munge_sumstats/*.sumstats.gz
	do
		if [ "`ls $GWASfile`" != "" ]
		then
			#echo $GWASfile
			trait=`echo "$GWASfile" | awk -F "/" '{print $NF }'`
			#echo ${trait}
			CGID=`echo "$trait" | awk -F "_" '{print $2 }'`
			#output=${output},${trait}
			echo '#!/bin/bash
#PBS -N '$CGID'_heritability
#PBS -q cu-1
#PBS -l walltime=1:00:00
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
#source activate ldsc 
module load R/3.6.2-anaconda3
source '$CUR_DIR'/.bashrc
source activate '$CUR_DIR'/.conda/envs/ldsc

~/software/ldsc/ldsc.py --chisq-max 80 --invert-anyway --return-silly-things --h2 '$GWASfile' --ref-ld-chr ~/software/ldsc/eur_w_ld_chr/ --w-ld-chr ~/software/ldsc/weights_hm3_no_hla/weights. --out ~/software/ldsc/result/20210729_ldsc_heritability/'$CGID'

wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

'> ~/software/ldsc/qsub/"$CGID".sh
					wait
					qsub ~/software/ldsc/qsub/"$CGID".sh
					wait	
	echo $CGID" is all finished"
	
	done
}

run_ldsc_genetic_correlation(){

mkdir -p $LDSC_output/ldsc_pair
	
for file1 in $LDSC_output/20210729_munge_sumstats/GLGC_*_result.txt.sumstats.gz
do
	i=`echo "$file1" | awk -F "/" '{print $NF }'| awk -F "_" '{print $2 }'`
	for file2 in $LDSC_output/20210729_munge_sumstats/GLGC_*_result.txt.sumstats.gz
	do
		j=`echo "$file2" | awk -F "/" '{print $NF }'| awk -F "_" '{print $2 }'`
		if [ "$i" != "$j" ] 
		then
			#echo $j
			cd $LDSC_output/20210729_ldsc_pair
			linearcorrName1=`echo "$i".vs."$j"_linearcorr.log`
			linearcorrName2=`echo "$j".vs."$i"_linearcorr.log`
			echo "$linearcorrName1"_"$linearcorrName2"
			if [[ ! -f "$linearcorrName1" ]] && [[ ! -f "$linearcorrName2" ]]
			then

				echo '#!/bin/bash
#PBS -N '$i'.vs.'$j'
#PBS -q cu-1
#PBS -l walltime=1:00:00
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
#source activate ldsc 
module load R/3.6.2-anaconda3
source '$CUR_DIR'/.bashrc
source activate '$CUR_DIR'/.conda/envs/ldsc
cd '$LDSC_output'/20210729_ldsc_pair
'$CUR_DIR'/software/ldsc/ldsc.py --chisq-max 80 --invert-anyway --return-silly-things --rg '$LDSC_output'/20210729_munge_sumstats/GLGC_'$i'_result.txt.sumstats.gz,'$LDSC_output'/20210729_munge_sumstats/GLGC_'$j'_result.txt.sumstats.gz --ref-ld-chr '$CUR_DIR'/software/ldsc/eur_w_ld_chr/ --w-ld-chr '$CUR_DIR'/software/ldsc/weights_hm3_no_hla/weights. --out '$LDSC_output'/20210729_ldsc_pair/'$i'.vs.'$j'_linearcorr
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

'> ~/software/ldsc/qsub/"$i".vs."$j".sh
					wait
					qsub ~/software/ldsc/qsub/"$i".vs."$j".sh
					wait
				
			echo $trait" is all finished"
			fi
		fi


	done
done
}

main
