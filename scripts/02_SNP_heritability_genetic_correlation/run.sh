main(){
      global_parameters
      creat_folder
      run_menge_sumstats
      run_ldsc_heritability
      run_ldsc_genetic_correlation
}

global_parameters(){
	#CUR_DIR=$(cd `dirname $0` && pwd -P)
	CUR_DIR="/lustre/home/ymhu"
	LDSC_GWASs=$CUR_DIR/software/ldsc/20210729_LDSC_input
	LDSC_output=$CUR_DIR/software/ldsc/result
}


create_folder(){
	mkdir -p input/
	mkdir -p $LDSC_GWASs
	mkdir -p $LDSC_output
	mkdir -p $CUR_DIR/software/ldsc/qusb
	mkdir -p $LDSC_output/20210729_ldsc_pair
}

run_munge_sumstats(){
	var=1
	mkdir -p $LDSC_output/20210729_munge_sumstats
	for file in {CG0051,CG0052,CG0102,CG0251,CG0255,CG0256,CG0258,CG0276,CG0361,CG0006,CG0007,CG0054,CG0083,CG0103,CG0104,CG0123,CG0124,CG0118,CG0133,CG0134,CG0154,CG0155,CG0259,CG0265,CG0266,CG0277,CG0055,CG0056,CG0116,CG0117,CG0057,CG0003,CG0312,CG0008,CG0012,CG0027,CG0030,CG0033,CG0042,CG0058,CG0068,CG0071,CG0084,CG0085,CG0089,CG0105,CG0125,CG0128,CG0119,CG0356,CG0308,CG0139,CG0142,CG0150,CG0156,CG0171,CG0172,CG0175,CG0178,CG0181,CG0184,CG0198,CG0201,CG0204,CG0359,CG0218,CG0237,CG0366,CG0260,CG0281,CG0291,CG0294,CG0301,CG0320,CG0323,CG0353,CG0038,CG0317,CG0362,CG0419,CG0011,CG0398,CG0388,CG0418,CG0077,CG0079,CG0080,CG0108-B,CG0191-B,CG0193-B,CG0381,CG0086,CG0087,CG0417,CG0435,CG0109,CG0127,CG0385,CG0130,CG0416,CG0421,CG0431,CG0432,CG0370,CG0371,CG0391,CG0357,CG0423,CG0136,CG0137,CG0386,CG0212,CG0239,CG0425,CG0436,CG0437,CG0146,CG0192,CG0194,CG0196,CG0197,CG0394,CG0397,CG0209,CG0422,CG0207,CG0208,CG0210,CG0211,CG0213,CG0214,CG0215,CG0424,CG0216,CG0360,CG0221,CG0222,CG0223,CG0224,CG0225,CG0226,CG0227,CG0228,CG0392,CG0393,CG0428,CG0429,CG0238,CG0387,CG0261,CG0268,CG0300,CG0304,CG0420,CG0327,CG0328,CG0342,CG0343,CG0344,CG0382,CG0383,CG0384,CG0426,CG0427,CG0434,CG0039,CG0345,CG0346,CG0364}
	do
		GWASfile=$LDSC_GWASs/GLGC_${file}_result.txt
	#for GWASfile in $LDSC_GWASs/GLGC_*
	#do
		#if [ "`ls $GWASfile`" != "" ]
		if [[ -f "$GWASfile" ]]
		then
			echo $GWASfile
			#echo $var
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
			#fi
			echo $trait" is all finished"
		fi
process=`qstat | grep "ymhu"|wc -l`
		while (( process >= 600))
		do
			echo "Current number of jobs is larger than 300"
			echo "Wait another 5 minutes"
			sleep 5m
			process=`qstat | grep "ymhu"|wc -l`
		done
		if [ $var -lt 1 ]
		then
			let "var+=1"
			continue
		fi
		if [ $var -gt 1000000 ]
		then
			break
		fi
		let "var += 1"
done

}
run_ldsc_heritability(){
	mkdir -p /lustre/home/ymhu/software/ldsc/result/20210729_ldsc_heritability
	for GWASfile in $LDSC_output/20210729_munge_sumstats/*.sumstats.gz
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
/lustre/home/ymhu/software/ldsc/ldsc.py --chisq-max 80 --invert-anyway --return-silly-things --h2 '$GWASfile' --ref-ld-chr /lustre/home/ymhu/software/ldsc/eur_w_ld_chr/ --w-ld-chr /lustre/home/ymhu/software/ldsc/weights_hm3_no_hla/weights. --out /lustre/home/ymhu/software/ldsc/result/20210729_ldsc_heritability/'$CGID'
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

'> /lustre/home/ymhu/software/ldsc/qsub/"$CGID".sh
					wait
					qsub /lustre/home/ymhu/software/ldsc/qsub/"$CGID".sh
					wait	
	echo $CGID" is all finished"
	fi
	process=`qstat | grep "ymhu"|wc -l`
		while (( process >= 600))
		do
			echo "Current number of jobs is larger than 300"
			echo "Wait another 1 minutes"
			sleep 1m
			process=`qstat | grep "ymhu"|wc -l`
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
}

run_ldsc_genetic_correlation(){
	var=1
	mkdir -p $LDSC_output/20210729_ldsc_pair
#for i in {CG0001,CG0002,CG0003,CG0004,CG0005,CG0006,CG0008,CG0009,CG0010,CG0011,CG0012,CG0013,CG0014,CG0015,CG0016,CG0017,CG0018,CG0019,CG0020,CG0021,CG0022,CG0023,CG0024,CG0025,CG0026,CG0027,CG0028,CG0029,CG0030,CG0031,CG0032,CG0033,CG0034,CG0035,CG0036,CG0037,CG0038,CG0039,CG0040,CG0041,CG0042,CG0043,CG0044,CG0045,CG0046,CG0048,CG0049,CG0051,CG0052,CG0054,CG0055,CG0056,CG0057,CG0058,CG0059,CG0060,CG0061,CG0062,CG0063,CG0064,CG0065,CG0066,CG0067,CG0068,CG0069,CG0070,CG0071,CG0072,CG0073,CG0074,CG0075,CG0076,CG0077,CG0078,CG0079,CG0080,CG0081,CG0083,CG0084,CG0085,CG0086,CG0087,CG0088,CG0089,CG0090,CG0091,CG0092,CG0093,CG0094,CG0095,CG0096,CG0097,CG0098,CG0099,CG0100,CG0103,CG0105,CG0106,CG0107,CG0108,CG0109,CG0110,CG0111,CG0112,CG0113,CG0114,CG0115,CG0116,CG0117,CG0118,CG0119,CG0120,CG0121,CG0123,CG0125,CG0126,CG0127,CG0128,CG0129,CG0130,CG0131,CG0133,CG0136,CG0137,CG0138,CG0139,CG0140,CG0141,CG0142,CG0143,CG0144,CG0145,CG0146,CG0147,CG0148,CG0149,CG0150,CG0154,CG0156,CG0157,CG0158,CG0159,CG0160,CG0161,CG0162,CG0163,CG0164,CG0165,CG0166,CG0167,CG0168,CG0169,CG0170,CG0171,CG0172,CG0173,CG0174,CG0175,CG0176,CG0177,CG0178,CG0179,CG0180,CG0181,CG0182,CG0183,CG0184,CG0185,CG0186,CG0187,CG0188,CG0189,CG0190,CG0191,CG0192,CG0193,CG0194,CG0195,CG0196,CG0197,CG0198,CG0199,CG0200,CG0201,CG0202,CG0203,CG0204,CG0205,CG0206,CG0207,CG0208,CG0209,CG0210,CG0211,CG0212,CG0213,CG0214,CG0215,CG0216,CG0217,CG0218,CG0219,CG0220,CG0221,CG0222,CG0223,CG0224,CG0225,CG0226,CG0227,CG0228,CG0229,CG0230,CG0231,CG0232,CG0233,CG0234,CG0235,CG0237,CG0238,CG0239,CG0245,CG0246,CG0247,CG0248,CG0249,CG0250,CG0251,CG0252,CG0253,CG0254,CG0255,CG0256,CG0257,CG0258,CG0259,CG0260,CG0261,CG0262,CG0265,CG0268,CG0269,CG0270,CG0271,CG0273,CG0276,CG0277,CG0281,CG0282,CG0283,CG0284,CG0285,CG0286,CG0287,CG0288,CG0289,CG0290,CG0291,CG0292,CG0293,CG0294,CG0295,CG0296,CG0297,CG0298,CG0299,CG0300,CG0301,CG0302,CG0303,CG0304,CG0305,CG0306,CG0307,CG0308,CG0309,CG0310,CG0311,CG0312,CG0313,CG0314,CG0315,CG0316,CG0317,CG0318,CG0319,CG0320,CG0321,CG0322,CG0323,CG0324,CG0325,CG0326,CG0327,CG0328,CG0329,CG0330,CG0331,CG0332,CG0333,CG0334,CG0335,CG0336,CG0337,CG0338,CG0339,CG0340,CG0341,CG0342,CG0343,CG0344,CG0345,CG0346,CG0347,CG0348,CG0349,CG0350,CG0351,CG0352,CG0353,CG0354,CG0355,CG0356,CG0357,CG0358,CG0359,CG0360,CG0362,CG0363,CG0364,CG0365,CG0366,CG0367,CG0368,CG0369,CG0370,CG0371,CG0372,CG0373,CG0374,CG0375,CG0376,CG0377,CG0378,CG0379,CG0380,CG0381,CG0382,CG0383,CG0384,CG0385,CG0386,CG0387,CG0388,CG0389,CG0390,CG0391,CG0392,CG0393,CG0394,CG0395,CG0396,CG0397}; do
#	for j in {CG0001,CG0002,CG0003,CG0004,CG0005,CG0006,CG0008,CG0009,CG0010,CG0011,CG0012,CG0013,CG0014,CG0015,CG0016,CG0017,CG0018,CG0019,CG0020,CG0021,CG0022,CG0023,CG0024,CG0025,CG0026,CG0027,CG0028,CG0029,CG0030,CG0031,CG0032,CG0033,CG0034,CG0035,CG0036,CG0037,CG0038,CG0039,CG0040,CG0041,CG0042,CG0043,CG0044,CG0045,CG0046,CG0048,CG0049,CG0051,CG0052,CG0054,CG0055,CG0056,CG0057,CG0058,CG0059,CG0060,CG0061,CG0062,CG0063,CG0064,CG0065,CG0066,CG0067,CG0068,CG0069,CG0070,CG0071,CG0072,CG0073,CG0074,CG0075,CG0076,CG0077,CG0078,CG0079,CG0080,CG0081,CG0083,CG0084,CG0085,CG0086,CG0087,CG0088,CG0089,CG0090,CG0091,CG0092,CG0093,CG0094,CG0095,CG0096,CG0097,CG0098,CG0099,CG0100,CG0103,CG0105,CG0106,CG0107,CG0108,CG0109,CG0110,CG0111,CG0112,CG0113,CG0114,CG0115,CG0116,CG0117,CG0118,CG0119,CG0120,CG0121,CG0123,CG0125,CG0126,CG0127,CG0128,CG0129,CG0130,CG0131,CG0133,CG0136,CG0137,CG0138,CG0139,CG0140,CG0141,CG0142,CG0143,CG0144,CG0145,CG0146,CG0147,CG0148,CG0149,CG0150,CG0154,CG0156,CG0157,CG0158,CG0159,CG0160,CG0161,CG0162,CG0163,CG0164,CG0165,CG0166,CG0167,CG0168,CG0169,CG0170,CG0171,CG0172,CG0173,CG0174,CG0175,CG0176,CG0177,CG0178,CG0179,CG0180,CG0181,CG0182,CG0183,CG0184,CG0185,CG0186,CG0187,CG0188,CG0189,CG0190,CG0191,CG0192,CG0193,CG0194,CG0195,CG0196,CG0197,CG0198,CG0199,CG0200,CG0201,CG0202,CG0203,CG0204,CG0205,CG0206,CG0207,CG0208,CG0209,CG0210,CG0211,CG0212,CG0213,CG0214,CG0215,CG0216,CG0217,CG0218,CG0219,CG0220,CG0221,CG0222,CG0223,CG0224,CG0225,CG0226,CG0227,CG0228,CG0229,CG0230,CG0231,CG0232,CG0233,CG0234,CG0235,CG0237,CG0238,CG0239,CG0245,CG0246,CG0247,CG0248,CG0249,CG0250,CG0251,CG0252,CG0253,CG0254,CG0255,CG0256,CG0257,CG0258,CG0259,CG0260,CG0261,CG0262,CG0265,CG0268,CG0269,CG0270,CG0271,CG0273,CG0276,CG0277,CG0281,CG0282,CG0283,CG0284,CG0285,CG0286,CG0287,CG0288,CG0289,CG0290,CG0291,CG0292,CG0293,CG0294,CG0295,CG0296,CG0297,CG0298,CG0299,CG0300,CG0301,CG0302,CG0303,CG0304,CG0305,CG0306,CG0307,CG0308,CG0309,CG0310,CG0311,CG0312,CG0313,CG0314,CG0315,CG0316,CG0317,CG0318,CG0319,CG0320,CG0321,CG0322,CG0323,CG0324,CG0325,CG0326,CG0327,CG0328,CG0329,CG0330,CG0331,CG0332,CG0333,CG0334,CG0335,CG0336,CG0337,CG0338,CG0339,CG0340,CG0341,CG0342,CG0343,CG0344,CG0345,CG0346,CG0347,CG0348,CG0349,CG0350,CG0351,CG0352,CG0353,CG0354,CG0355,CG0356,CG0357,CG0358,CG0359,CG0360,CG0362,CG0363,CG0364,CG0365,CG0366,CG0367,CG0368,CG0369,CG0370,CG0371,CG0372,CG0373,CG0374,CG0375,CG0376,CG0377,CG0378,CG0379,CG0380,CG0381,CG0382,CG0383,CG0384,CG0385,CG0386,CG0387,CG0388,CG0389,CG0390,CG0391,CG0392,CG0393,CG0394,CG0395,CG0396,CG0397}; do
for file1 in $LDSC_output/20210729_munge_sumstats/GLGC_*_result.txt.sumstats.gz
do
	i=`echo "$file1" | awk -F "/" '{print $NF }'| awk -F "_" '{print $2 }'`
	for file2 in $LDSC_output/20210729_munge_sumstats/GLGC_*_result.txt.sumstats.gz
	do
		j=`echo "$file2" | awk -F "/" '{print $NF }'| awk -F "_" '{print $2 }'`
		if [ "$i" != "$j" ] #不等于-ne
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

'> /lustre/home/ymhu/software/ldsc/qsub/"$i".vs."$j".sh
					wait
					qsub /lustre/home/ymhu/software/ldsc/qsub/"$i".vs."$j".sh
					wait
					#fi
				#fi
			#fi
			echo $trait" is all finished"
			fi
		fi

	process=`qstat | grep "ymhu"|wc -l`
		while (( process >= 300))
		do
			echo "Current number of jobs is larger than 300"
			echo "Wait another 1 minutes"
			sleep 1m
			process=`qstat | grep "ymhu"|wc -l`
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
done
}

main
