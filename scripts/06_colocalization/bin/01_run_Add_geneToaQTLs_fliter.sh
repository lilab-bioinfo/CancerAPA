#!/bin/bash
#PBS -N lilei_lab
#PBS -q fat-1
#PBS -l nodes=1:ppn=2
#PBS -V
#PBS -S /bin/bash
#cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

cd /lustre/home/ymhu/aQTL_coloc
module load R/3.6.2-anaconda3
mkdir -p /lustre/home/ymhu/aQTL_coloc/input/20210728_aQTLs_postionTohg19_for_coloc
python /lustre/home/ymhu/aQTL_coloc/bin/Add_geneToaQTLs_fliter.py -b /lustre/home/ymhu/aQTL_coloc/bin/hg38_Refseq_id_from_UCSC.txt -p /lustre/home/ymhu/aQTL_coloc/input/20210718_aQTLs_postionTohg19_for_wenyan -o /lustre/home/ymhu/aQTL_coloc/input/20210728_aQTLs_postionTohg19_for_coloc
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/


