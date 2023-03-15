
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
#module load R/3.6.2-anaconda3
mkdir -p /lustre/home/ymhu/aQTL_coloc/output/aQTLs
for file in /lustre/home/ymhu/aQTL_coloc/input/20210728_aQTLs_postionTohg19_for_coloc/*
	do
		tissueName=`echo "$file" |awk -F"/" '{print $NF;exit}'|awk -F".cis_aqtl_hg19" '{print $1;exit}'`
		cd /lustre/home/ymhu/aQTL_coloc/output/aQTLs
		if [[ ! -d "$tissueName" ]]
			then
cd /lustre/home/ymhu/aQTL_coloc/input/20210728_aQTLs_postionTohg19_for_coloc/
mkdir -p /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}
head -n1 $file >/lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_temptitle.txt
tail -n +2 $file |sort -k3 > /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_tempsorted.txt
wait
cat /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_temptitle.txt /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_tempsorted.txt > /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_tempmerged.txt
wait
python2 /lustre/home/ymhu/coloc_analysis/bin/extract_gene.py /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_tempmerged.txt /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}
wait
rm /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_temptitle.txt
rm /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_tempsorted.txt
rm /lustre/home/ymhu/aQTL_coloc/output/aQTLs/${tissueName}_tempmerged.txt

wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
#module unload R/3.6.2-anaconda3
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait

fi
echo $tissueName" is finished"
done
