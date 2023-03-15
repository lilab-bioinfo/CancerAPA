#!/bin/bash
#PBS -N Brain_Frontal_Cortex_BA9
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
cd /lustre/home/ymhu/aQTL_pipeline
python /lustre/home/ymhu/aQTL_pipeline/bin/aQTL_format_PAINTOR.py -c /lustre/home/ymhu/aQTL_pipeline/output_QTL/Brain_Frontal_Cortex_BA9.cis_eqtl_all_approachb.txt.gz -i /lustre/home/ymhu/TWAS/input/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz -o /lustre/home/ymhu/aQTL_pipeline
python /lustre/home/ymhu/aQTL_coloc/bin/Add_geneToaQTLs_fliter.py -b /lustre/home/ymhu/aQTL_coloc/bin/hg38_Refseq_id_from_UCSC.txt -p /lustre/home/ymhu/aQTL_coloc/input/20210718_aQTLs_postionTohg19_for_wenyan -o /lustre/home/ymhu/aQTL_coloc/input/20210728_aQTLs_postionTohg19_for_coloc
mkdir -p /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9
cd /lustre/home/ymhu/aQTL_coloc/input/20210728_aQTLs_postionTohg19_for_coloc/
head -n1 Brain_Frontal_Cortex_BA9.cis_aqtl_hg19.txt >/lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_temptitle.txt
tail -n +2 Brain_Frontal_Cortex_BA9.cis_aqtl_hg19.txt |sort -k3 > /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_tempsorted.txt
wait
cat /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_temptitle.txt /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_tempsorted.txt > /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_tempmerged.txt
wait
python2 /lustre/home/ymhu/coloc_analysis/bin/extract_gene.py /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_tempmerged.txt /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9
wait
rm /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_temptitle.txt
rm /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_tempsorted.txt
rm /lustre/home/ymhu/aQTL_coloc/output/aQTLs/Brain_Frontal_Cortex_BA9_tempmerged.txt
wait
echo "processs will sleep 5s"
sleep 5
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
