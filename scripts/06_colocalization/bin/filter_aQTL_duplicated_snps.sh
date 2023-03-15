#!/bin/bash
#PBS -N aqtl_filter
#PBS -q cu-1
#PBS -l nodes=1:ppn=2
#PBS -V
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$
#gene    chr     position        rsid    pval_nominal    beta    varbeta N       MAF
#chrposition	gene	p-value	N	rsid
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
module load R/3.6.2-anaconda3
for file in /lustre/home/ymhu/aQTL_coloc/output/aQTLs_20210728/*
do
tissueName=`echo "$file" |awk -F"/" '{print $NF;exit}'`
for gene in /lustre/home/ymhu/aQTL_coloc/output/aQTLs_20210728/"$tissueName"/*.txt
do
geneName=`echo "$gene" |awk -F"/" '{print $NF;exit}'`
python /lustre/home/ymhu/sQTL/bin/filter_QTL_duplicated_snps.py -i /lustre/home/ymhu/aQTL_coloc/output/aQTLs_20210728/"$tissueName"/"$geneName" -o /lustre/home/ymhu/aQTL_coloc/output/aQTLs_20210728/"$tissueName"/"$geneName"
done
echo $tissueName" is finished"
done
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
module unload R/3.6.2-anaconda3
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/