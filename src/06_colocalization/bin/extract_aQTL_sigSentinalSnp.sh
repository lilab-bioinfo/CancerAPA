#!/bin/bash
#PBS -N aQTL_sigSentinalSnp
#PBS -q fat-1
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
cd ~/aQTL_coloc/bin

module load R/3.6.2-anaconda3

python sigSentinalSnp_all.py
wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

module unload R/3.6.2-anaconda3
wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
