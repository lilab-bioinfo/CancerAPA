main(){
        global_parameters
        get_tisssue_gene_list
        split_file
        finemapping_analysis
}

global_parameters(){
        DIR=`pwd`
}

get_tissue_gene_list(){

      echo '#!/bin/bash
#SBATCH --job-name=\'get_list\'
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --error=%j.err
#SBATCH --output=%j.out
        

        for file in ~/Data/GTEx_v8/aQTL/hg38/*.cis_aqtl_all_approachb.txt.gz
        do
        filename=${file##*/}
        tissue=${filename%.cis_aqtl_all_approachb.txt.gz}
        echo ${tissue}
        less ${file} | awk \'$NF<0.05\' | cut -f2 |sort|uniq|while read i; do echo -e "${tissue}\t${i}";done > $DIR/input/${tissue}_gene_list.txt
        done' > get_tissue_gene_list.sh
wait    

sbatch get_tissue_gene_list.sh  

}

split_file(){
for file in ~/Data/GTEx_v8/aQTL/hg38/*.cis_aqtl_all_approachb.txt.gz
                do
                filename=${file##*/}
                tissue=${filename%.cis_aqtl_all_approachb.txt.gz}
                echo ${tissue}
                echo -e '#!/bin/bash
#SBATCH --job-name=\'get_list\'
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --error=%j.err
#SBATCH --output=%j.out

echo "start at:"
date

~/miniconda2/bin/python ${DIR}/bin/Split_file.caviar.py '${tissue}'

echo "ends at:"
date

' > $DIR/split/${tissue}_r.sh

}

finemapping_analysis(){
VCF="~/Data/GTEx_v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz" 
GENELOC="~/aQTL_pipeline/Input_gene_locs/gene_3UTR_loc_hg38.txt" ## hg19

for LSB_JOBINDEX in `seq 1 189116` 
        do
                EGENE=`less $DIR/input/GTEx_gene_list.txt|sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p|cut -f 2 `
                TISSUE=`less $DIR/input/GTEx_gene_list.txt|sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p|cut -f 1 `
                echo ${LSB_JOBINDEX}":"$EGENE":"$TISSUE
                mkdir -p $DIR/output/${TISSUE}/${EGENE}/caviar/
                echo "#!/bin/bash
#SBATCH --job-name=\'caviar\'
#SBATCH --partition=fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --error=%j.err
#SBATCH --output=%j.out

##################################

echo "start at:"
date

module load java/jdk-1.8.0_241
module load anaconda/anaconda2-2019.10-py27

cd $DIR

EGENE=\"`less $DIR/input/GTEx_gene_list.txt|sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p|cut -f 2 `\"
TISSUE=\"`less $DIR/input/GTEx_gene_list.txt|sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p|cut -f 1 `\"
REGION=\`cat ${GENELOC} | awk -v EGENE=\${EGENE} -v SLOP=1000000 '{ if (\$1==EGENE) { POS_START=\$3-SLOP; POS_END=\$4+SLOP; if (POS_START<0) POS_START=0; print \$2\":\"POS_START\"-\"POS_END } }' \`

mkdir -p $DIR/output/\${TISSUE}/\${EGENE}/mott
tabix -h $VCF \$REGION >  $DIR/output/\${TISSUE}/\${EGENE}/vcf.tmp

mkdir -p $DIR/output/\${TISSUE}/\${EGENE}/caviar
python $DIR/bin/ld_continuous.py \
-i $DIR/output/\${TISSUE}/\${EGENE}/vcf.tmp \
-v $DIR/output/\${TISSUE}/\${EGENE}/\${TISSUE}_3aQTL_all.snp \
-f GT \
-a r \
> $DIR/output/\${TISSUE}/\${EGENE}/caviar/\${TISSUE}_3aQTL_ld

/lustre/software/caviar/caviar/CAVIAR-C++/CAVIAR \
-o $DIR/output/\${TISSUE}/\${EGENE}/caviar/\${TISSUE}_3aQTL.cts \
-l $DIR/output/\${TISSUE}/\${EGENE}/caviar/\${TISSUE}_3aQTL_ld \
-z $DIR/output/\${TISSUE}/\${EGENE}/\${TISSUE}_3aQTL_all_ttest.txt \
-r 0.95 \
-c 2 \
-f 1

echo "end at:"
date
        " > $DIR/caviar_analysis/job_${LSB_JOBINDEX}.sh
        wait
        
        cd $DIR/output/${TISSUE}/${EGENE}/caviar
        
        if [[ -s "${TISSUE}_3aQTL.cts_set" ]];
        then
                echo "file exists"
                continue
        else
                echo "analysing"
                cd $DIR/caviar_analysis/log
                sbatch $DIR/caviar_analysis/job_${LSB_JOBINDEX}.sh &
                wait
                
done
}

main
