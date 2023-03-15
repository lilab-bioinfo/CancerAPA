main(){ 
  global_parameters
  compute_weights
  generate_summary
  run_3aTWAS_analysis
  conditional_analysis
}

global_parameters(){
    WORK_DIR=`pwd`
    INPUT_DIR=$WORK_DIR/input
    OUTPUT_DIR=$WORK_DIR/output
}

compute_weights(){
for tiss in `less tissue.list`
  do

    PRE="${tiss}"
    
    mkdir -p $WORK_DIR/slurm/${tiss}
    
    PRE_GECP="${PRE}_combined_All_PDUIs_clean.txt"
    BATCH_START=1
    BATCH_END=`less ${INPUT_DIR}/APA_matrix/$PRE_GERP}|wc -l`
    NR="${BATCH_START}_${BATCH_END}"
    
    less $INPUT_DIR/APA_matrix/$PRE_GEXP |awk -vs=$BATCH_START -ve=$BATCH_END 'NR > s && NR <= e' |  while read PARAM; do
    GNAME=`echo $PARAM | awk '{ print $4 }' |awk -F "|" '{print $1}'`
    gene=`echo $PARAM | awk '{ print $4 }' |awk -F "|" '{print $2}'`
    chr=`echo $PARAM | awk '{ print $4 }' |awk -F "|" '{print $3}'`
    strand=`echo $PARAM | awk '{ print $4 }' |awk -F "|" '{print $4}'`
    
    Name=`echo ${GNAME}"@"${gene}"@"${chr}"@"${strand}`
    CHR=`echo $PARAM | awk '{ print $1 }' |awk -F "chr" '{print $2}'`
    P0=`echo $PARAM | awk '{ p=$2 - 500e3; if(p<0) p=0; print p; }'`
    P1=`echo $PARAM | awk '{ print $3 + 500e3 }'`
    
    WEIGHT="$OUTPUT_DIR/WEIGHTS/${PRE}/${PRE}.${Name}.wgt.RDat"
    
    echo '#!/bin/bash
#SBATCH --job-name='"'"''$tiss'_'$GNAME'_build_3a_weight'"'"'
#SBATCH --partition=cpuPartition
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --error=%j.err
#SBATCH --output=%j.out

##################################

echo "process will start at : "

date

echo "++++++++++++++++++++++++++++++++++++++++"

GCTA="~/biosoft/gcta_1.93.3beta2/gcta64"
PLINK="~/biosoft/plink_linux_x86_64_20210606/plink"
GEMMA="~/biosoft/gemma-0.98.5-linux-static-AMD64"
FUSION="~/fusion/bin/fusion_twas-master"

PRE="'${tiss}'"

WORK_DIR="pwd"
INPUT_DIR=$WORK_DIR/input
OUTPUT_DIR=$WORK_DIR/output

cd $OUTPUT_DIR

### Input files needed:
# "$PRE.v7.normalized_expression.bed.gz.HEADER" is generated by `cat [bed.gz] | head -n1 | tr '\t' '\n'`
# "$PRE.v7.covariates.txt.covar" is a PLINK format covariates file
# "$PRE.v7.normalized_expression.bed.gz" is the gzipped bed format expression matrix

## PATH to GTEx v8 APA MATRIX:
PRE_GEXP="${PRE}_combined_All_PDUIs_clean.txt"

less $INPUT_DIR/APA_matrix/$PRE_GEXP | head -n1 | tr '"'"' '"'"' '"'"'\n'"'"' | tail -n+5 | awk '"'"'{ print $1,$1 }'"'"' >$OUTPUT_DIR/HEADER/${PRE_GEXP}.HEADER

BATCH_START=1
BATCH_END=`less ${INPUT_DIR}/APA_matrix/${PRE_GEXP}|wc -l`

NR="${BATCH_START}_${BATCH_END}"
echo $NR

cd $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/tmp/${PRE}.${NR}
mkdir -p $OUTPUT_DIR/HSQ

COVAR="$INPUT_DIR/input/covariate_plink_format/${PRE}.pdui.peer.covariates.txt.bak"

OUT="$OUTPUT_DIR/tmp/${PRE}.'${NR}'/${PRE}.'${GNAME}'"

echo '$GNAME' '$CHR' '$P0' '$P1'
echo "'$PARAM'" |tr '"'"' '"'"' '"'"'\n'"'"' |tail -n+5> $OUTPUT_DIR/tmp/${PRE}.'${NR}'/${PRE}.'${GNAME}'.temp.pheno
paste $OUTPUT_DIR/HEADER/${PRE_GEXP}.HEADER $OUTPUT_DIR/tmp/${PRE}.'${NR}'/${PRE}.'${GNAME}'.temp.pheno > ${OUT}.pheno

rm $OUTPUT_DIR/tmp/${PRE}.${NR}/${PRE}.'${GNAME}'.temp.pheno

# Extract the locus around this gene
$PLINK --allow-no-sex --silent --bfile ${WORK_DIR}/genotypes.bed.rsid/GTEx_v8.${PRE} --chr '$CHR' --from-bp '$P0' --to-bp '$P1' --make-bed --out $OUT --pheno ${OUT}.pheno --keep $OUTPUT_DIR/HEADER/${PRE_GEXP}.HEADER --maf 0.01

# Run FUSION
module load R/3.6.2-anaconda3
mkdir -p $OUTPUT_DIR/WEIGHTS
mkdir -p $OUTPUT_DIR/WEIGHTS/${PRE}
FINAL_OUT="$OUTPUT_DIR/WEIGHTS/${PRE}/${PRE}.'${GNAME}'"

Rscript $FUSION/FUSION.compute_weights.R \
--bfile ${OUT} --tmp ${OUT}.tmp --out $FINAL_OUT --verbose 2 --save_hsq --PATH_plink  $PLINK --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,enet --covar ${COVAR} --hsq_p 0.05

cat ${FINAL_OUT}.hsq | awk -vw="${PRE} $w" '"'"'{ print w,$0 }'"'"' >> $OUTPUT_DIR/HSQ/${PRE}.${NR}.hsq

rm -f $FINAL_OUT.hsq
rm ${OUT}.*
    
wait
echo "process end at : "
date
    
wait' > $WORK_DIR/slurm/${tiss}/${GNAME}_compute_weight.sh

        sbatch $WORK_DIR/slurm/${tiss}/${GNAME}_compute_weight.sh
done
}

generate_summary(){
cd $OUTPUT_DIR/WEIGHTS
ls|while read line;do ls ${line}/* >`echo ${line}`.list; done

module load R/3.6.2-anaconda3

makeProfile(){
      file=$1
      tiss=`echo $file|awk -F "/" '{print $NF}'|awk -F "." '{print $1}'`
      cd $OUTPUT_DIR/WEIGHTS
      Rscript $fusion/bin/fusion_twas-master/utils/FUSION.profile_wgt.R $OUTPUT_DIR/WEIGHTS/${tiss}.list > $OUTPUT_DIR/WEIGHTS/${tiss}.profile 2 > $OUTPUT_DIR/WEIGHTS/${tiss}.profile.err
}

tissALL=(`ls $OUTPUT_DIR/WEIGHTS/*.list`)
for(( t=0; t<49; t++)){ makeProfile "${tissALL[$t]}" &}

}

run_3aTWAS_analysis(){

for file in  `ls $WORK_DIR/GWAS/*.sumstats`
      do
      GWAS=`echo $file|awk -F "/" '{print $NF}' |awk -F "_result" '{print $1}'`
      for tissueName in `less $WORK_DIR/bin/tissue.list`
           do
           mkdir -p $WORK_DIR/results/${GWAS}/${tissueName}
           echo '#!/bin/bash
#SBATCH --job-name=\'3aTWAS_'${GWAS}'\'
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=%j.err
#SBATCH --output=%j.out
      
##################################

module load R/3.6.2-anaconda3
for chr  in `less $WORK_DIR/WEIGHTS/'${tissueName}'.pos|awk \'NR>1{print $3}\'|sort -u`
                        do
                        outfile="$WORK_DIR/results/'${GWAS}'/'${tissueName}'/'${GWAS}'.'${tissueName}'.chr${chr}.dat"
                if [ -f ${outfile} ];then
                        echo "'${GWAS}'_"'${tissueName}'"_chr"${chr}" have done!"
                        continue
                else
                Rscript /lustre/home/hchen/2021-10-31-cancer-GWAS/aTWAS/fusion/bin/fusion_twas-master/FUSION.assoc_test.R \
                --sumstats '${file}' \
                --weights /lustre/home/hchen/2021-10-31-cancer-GWAS/aTWAS/fusion/2022-02-17-build-weight/output_FDR0.05_new/WEIGHTS/'${tissueName}'.pos  \
                --weights_dir /lustre/home/hchen/2021-10-31-cancer-GWAS/aTWAS/fusion/2022-02-17-build-weight/output_FDR0.05_new/WEIGHTS/ \
                --ref_ld_chr /lustre/home/hchen/2021-10-31-cancer-GWAS/aTWAS/fusion/bin/fusion_twas-master/LDREF/1000G.EUR. \
                --chr ${chr} \
                --out /lustre/home/hchen/2021-10-31-cancer-GWAS/aTWAS/fusion/results_new/'${GWAS}'/'${tissueName}'/'${GWAS}'.'${tissueName}'.chr${chr}.dat 
                
                fi
               
                done' > $WORK_DIR/slurm/${GWAS}/${GWAS}_${tissueName}_run_aTWAS.sh
                
                mkdir -p $WORK_DIR/slurm/${GWAS}/log
                cd $WORK_DIR/slurm/${GWAS}/log
                sbatch $WORK_DIR/slurm/${GWAS}/${GWAS}_${tissueName}_run_aTWAS.sh
                done
      done


}

contidional_analysis(){
###### filter the top results ###########
for tissueName in `less $WORK_DIR/bin/tissue.list`
    do
    num_weights=`less $WORK_DIR/WEIGHTS/"${tissueName}".list|wc -l`
    cd $WORK_DIR/results
    for GWAS in `ls $WORK_DIR/results/`
        do
        ls $WORK_DIR/results/${GWAS}/${tissueName}/*.dat|while read line;do cat $line |awk 'NR==1 || $NF <0.05/'${num_weights}'' >  `echo ${line}|awk -F ".dat" '{print $1}'`.top;done
    done
done

########## post_process #############
echo '#!/bin/bash
#SBATCH --job-name=\'post_process\'
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=%j.err
#SBATCH --output=%j.out
      
##################################

module load R/3.6.2-anaconda3
for file in `ls $WORK_DIR/GWASs/*.sumstats`
    do
    GWAS=`echo $file|basename |awk -F "_result" '{print $1}'`
        for CHR in `less $OUTPUT_DIR/WEIGHTS/${tissueName}.pos|awk \'NR>1{print $3}\'|sort -u`
        do
        Rscript $fusion/bin/fusion_twas-master/FUSION.post_process.R --sumstats ${file} \
        --input $WORK_DIR/resulsts/${GWAS}/${tissueName}/${GWAS}.${tissueName}.chr${CHR}.top \
        --out /lustre/home/hchen/2021-10-31-cancer-GWAS/aTWAS/fusion/results/${GWAS}/${tissueName}/${GWAS}.${tissueName}.chr${CHR}.top.analysis \
        --ref_ld_chr $fusion/bin/fusion_twas-master/LDREF/1000G.EUR \
        --chr $CHR \
        --plot \
        --locus_win 100000    
        done
done


' > $WORK_DIR/slurm/post_process.sh
    sbatch $WORK_DIR/slurm/post_process.sh
}


main
