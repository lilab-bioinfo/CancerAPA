main(){
	set_parameters
	create_folder
	#run_filterVCF
	#run_pastefilterVCF
	#run_filterVCF2value
	#run_peer
	run_matrixeqtl
	#boxplot_3aQTL
}

set_parameters(){
	tag='Skin_Not_Sun_Exposed_Suprapubic'
}

create_folder(){
	mkdir -p input output bin

}
run_filterVCF(){
	var=1
	mkdir -p ~/aQTL_pipeline/qsub/filterVCF
	mkdir -p ~/aQTL_pipeline/input/vcf
	mkdir -p ~/aQTL_pipeline/input/vcf/temp
	for i in {2..839}
	do
	cd ~/aQTL_pipeline/input/vcf/temp

	echo "#!/bin/bash
#PBS -N filterVCF_'$i'
#PBS -q cu-1
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -S /bin/bash

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
module load R/3.6.2-anaconda3
cd ~/aQTL_pipeline/input
#zcat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP.vcf.gz | cut -f'$i' - > ~/aQTL_pipeline/input/vcf/temp/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP_'$i' 
python3 ~/aQTL_pipeline/bin/GTEx_V8_single_filter_SNPs_addcovariates.py -p ~/aQTL_pipeline/input/vcf/temp/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP_'$i' -i ~/aQTL_pipeline/input/GTEx_covariates.v8.txt -o ~/aQTL_pipeline/input/filter/
wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

" > ~/aQTL_pipeline/qsub/filterVCF/filterVCF_${i}.sh
wait
qsub ~/aQTL_pipeline/qsub/filterVCF/filterVCF_${i}.sh
echo "filterVCF_"${i}" is all finished"

wait


done

}

run_pastefilterVCF(){
	mkdir -p ~/aQTL_pipeline/qsub/pastefilterVCF
	mkdir -p ~/aQTL_pipeline/output
	for i in ~/aQTL_pipeline/Input_Covariates/Covariates_by_tissue/*.v8.covariates.txt
	do
	traitName=`echo "$i" |awk -F"/" '{print $NF;exit}'|awk -F".v8." '{print $1;exit}'`
	cd ~/aQTL_pipeline/output
	if [[ ! -f "$traitName"_SNPs.vcf_genotype.vcf ]]
	then
	echo "#!/bin/bash
#PBS -N paste_'$traitName'
#PBS -q cu-1
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -S /bin/bash

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

module load R/3.6.2-anaconda3

cd ~/aQTL_pipeline/input/filter/

paste -d'	' id.txt ${traitName}_*.txt > ~/aQTL_pipeline/output/${traitName}_SNPs.vcf_genotype.vcf

wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

" > ~/aQTL_pipeline/qsub/pastefilterVCF/pastefilterVCF_${traitName}.sh

wait
qsub ~/aQTL_pipeline/qsub/pastefilterVCF/pastefilterVCF_${traitName}.sh
echo "filterVCF_"${traitName}" is all finished"

fi
wait

done


}
run_filterVCF2value(){
	mkdir -p ~/aQTL_pipeline/qsub/filterVCF
	mkdir -p ~/aQTL_pipeline/output_filter
	for i in ~/aQTL_pipeline/output/*_SNPs.vcf_genotype.vcf
	do
	traitName=`echo "$i" |awk -F"/" '{print $NF;exit}'`
	cd ~/aQTL_pipeline/output_filter
	if [[ ! -f "$traitName"_SNPs.vcf_genotype.vcf ]]
	then
	echo "#!/bin/bash
#PBS -N '$traitName'
#PBS -q cu-1
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -S /bin/bash

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
module load R/3.6.2-anaconda3
cd ~/aQTL_pipeline/output_filter
python3 ~/aQTL_pipeline/bin/GTEx_V8_filter_VCF_Genetype2value.py -p ${i} -o ~/aQTL_pipeline/output_filter/${traitName}
wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

" > ~/aQTL_pipeline/qsub/filterVCF/filterVCF_${traitName}.sh

wait
qsub ~/aQTL_pipeline/qsub/filterVCF/filterVCF_${traitName}.sh
echo "filterVCF_"${traitName}" is all finished"
fi
wait


done


}

run_peer(){
#Rscript bin/3UTR_QTL_peer.R ${tissue}" > ${tissue}_peer.sh

	for file in input/${tag}*combined_All_PDUIs_clean.txt
	do
	filename=${file##*/}
	tissue=${filename%_combined_All_PDUIs_clean.txt}
	cd ~/aQTL_pipeline/peer_pdui
	if [[ ! -f "$tissue".pdui.impute.txt.v2 ]]
	then
	echo "#!/bin/bash
#PBS -N '$filename'_peer
#PBS -q cu-1
#PBS -l nodes=1:ppn=1
#PBS -V
#PBS -S /bin/bash

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
module load R/3.6.2-anaconda3
Rscript ~/aQTL_pipeline/bin/3UTR_QTL_peer_impute.R ${tissue}
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

" > ~/aQTL_pipeline/qsub/${tissue}_peer.sh
wait
qsub ~/aQTL_pipeline/qsub/${tissue}_peer.sh
fi
wait
	done

#ls *expr.peer.txt|while read i; do echo -e "$i\t`awk '{print NF}' $i|head -n1`\t`awk '{print NF}' ../input/${i%.expr.peer.txt}_combined_All_PDUIs_clean.txt|head -n1`";done|less
}

run_matrixeqtl(){
#Rscript bin/run_matrixEQTL_without_plot_control_gene_exprs.R ${tissue}"> ${tissue}_v3.sh
	mkdir -p ~/aQTL_pipeline/output_QTL
	for file in output_filter/${tag}*_SNPs.vcf_genotype.vcf
	do
	filename=${file##*/}
	tissue=${filename%_SNPs.vcf_genotype.vcf}
	cd ~/aQTL_pipeline
	output=output_QTL/${tissue}.cis_eqtl_all_approachb.txt
	if [ ! -f ${output} ]
	then
	echo ${tissue}
	echo "#!/bin/bash
#PBS -N '$tissue'_matrixeqtl
#PBS -q fat-1
#PBS -l nodes=1:ppn=4
#PBS -V
#PBS -S /bin/bash

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
module load R/3.6.2-anaconda3
cd ~/aQTL_pipeline
Rscript bin/run_matrixEQTL_without_plot_covariates_final.R ${tissue}
wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
"> ~/aQTL_pipeline/qsub/${tissue}_v3.sh

wait
qsub ~/aQTL_pipeline/qsub/${tissue}_v3.sh
wait

	fi
	done
}

boxplot_3aQTL(){
	module load R/3.4.1

	#Figure S7
	#for file in output_QTL/*all.txt
	#do
	#filename=${file##*/}
	#tissue=${filename%.cis_eqtl_all.txt}
	#Rscript bin/QTL_plot.R -t ${tissue} -s chr7_128589427_G_A -g "NM_001098630|IRF5|chr7|+"
	#done

	#Rscript bin/QTL_plot.R -t Muscle_Skeletal -s chr16_2088312_C_T -g "NM_004785|SLC9A3R2|chr16|+"
	#Rscript bin/QTL_plot.R -t Muscle_Skeletal -s chr2_69550036_T_C -g "NM_002056|GFPT1|chr2|-"
	#Rscript bin/QTL_plot.R -t Thyroid -s chr19_18685964_G_T -g "NM_003333|UBA52|chr19|+"
	#Rscript bin/QTL_plot.R -t Muscle_Skeletal -s chr16_84156280_C_G -g "NM_001146051|HSDL1|chr16|-"
	#Rscript bin/QTL_plot.R -t Whole_Blood -s chr1_40219065_C_T -g "NM_006112|PPIE|chr1|+"

	#rs66534072
	Rscript bin/QTL_plot.R -t Whole_Blood -s chr22_21936152_C_G -g "NR_046082|UBE2L3|chr22|+"

	#
	Rscript bin/QTL_plot.R -t Liver -s chr12_109925489_A_G -g "NM_052845|MMAB|chr12|-"

	#
	Rscript bin/QTL_plot.R -t Muscle_Skeletal -s chr16_4015046_G_A -g "NM_001116|ADCY9|chr16|-"

	#
	Rscript bin/QTL_plot.R -t Cells_Transformed_fibroblasts -s chr5_131827775_A_G -g "NM_002198|IRF1|chr5|-"


	#Rscript bin/QTL_plot.R -t Muscle_Skeletal -s chr16_4015046_G_A -g "NM_001116|ADCY9|chr16|-"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr1_119555759_C_T -g "NM_201263|WARS2|chr1|-"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr1_151839838_A_C -g "NM_053055|THEM4|chr1|-"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr1_160968005_G_C -g "NM_016946|F11R|chr1|-"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr3_150302010_C_T -g "NM_032025|EIF2A|chr3|+"

	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr12_51138862_C_T -g "NM_173602|DIP2B|chr12|+"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr7_128589427_G_A -g "NM_001098630|IRF5|chr7|+"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr7_128630313_T_C -g "NM_001098630|IRF5|chr7|+"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr6_32611100_G_A -g "NM_002122|HLA-DQA1|chr6|+"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr10_16555528_A_T -g "NM_030664|PTER|chr10|+"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr22_22048862_C_T -g "NM_014337|PPIL2|chr22|+"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr5_96111371_G_A -g "NM_001040458|ERAP1|chr5|-"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr5_96082433_T_C -g "NM_001040458|ERAP1|chr5|-"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr5_112359020_A_C -g "NM_001242377|DCP2|chr5|+"
	#Rscript bin/QTL_plot.R -t Cells_EBV-transformed_lymphocytes -s chr16_56393191_T_C -g "NM_001144|AMFR|chr16|-"


	#Rscript bin/QTL_plot.R -t Muscle_Skeletal -s chr16_2088312_C_T -g "NM_004785|SLC9A3R2|chr16|+"
	#Rscript bin/QTL_plot.R -t Muscle_Skeletal -s chr2_69550036_T_C -g "NM_002056|GFPT1|chr2|-"
	#Rscript bin/QTL_plot.R -t Thyroid -s chr19_18685964_G_T -g "NM_003333|UBA52|chr19|+"
	#Rscript bin/QTL_plot.R -t Muscle_Skeletal -s chr16_84156280_C_G -g "NM_001146051|HSDL1|chr16|-"
	#Rscript bin/QTL_plot.R -t Whole_Blood -s chr1_40219065_C_T -g "NM_006112|PPIE|chr1|+"

	#Rscript bin/QTL_plot.R -t Whole_Blood -s chr2_70488470_C_T -g "NM_017567|NAGK|chr2|+"


}


main

