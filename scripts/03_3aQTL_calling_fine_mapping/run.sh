main(){
	set_parameters
	create_folder
	#run_filterVCF
	#run_pastefilterVCF
	#run_filterVCF2value
	#link_all_pdui
	#get_id_file
	#filter_vcftools
	#get_BioAlcidae
	#vcf_to_snpmatrix
	#run_peer
	run_matrixeqtl
	#convert_to_torus
	#convert_to_torus_2
	#compare_with_other_covariate
	#compare_with_other_covariate_2
	#boxplot_3aQTL
	#generate_qtl_plots
	#generate_trans_qtl_plots
	#get_sig_qtl_list
	#generate_subset_genotype
	#run_permutation_test
	#get_snp_count
	#get_final_qtl
	#get_final_qtl_v2
	#eGene_filter
	#run_vcf_filtering
	#add_rsid
	#plot_specific_genes
	#find_lead_3aQTL
	#find_lead_3aQTL_2
	#compare_best_variants
	#convert_3QTL_to_bed
	#3QTLs_to_publish
	#compare_qtl_differenes
}

set_parameters(){
	#Heart_Left_Ventricle
	#Skin_Not_Sun_Exposed_Suprapubic
	#Artery_Tibial
	#tag='Artery_Tibial'
	tag=''

}

create_folder(){
	mkdir -p input output bin

}
run_filterVCF(){
	var=1
	mkdir -p /lustre/home/ymhu/aQTL_pipeline/qsub/filterVCF
	mkdir -p /lustre/home/ymhu/aQTL_pipeline/input/vcf
	mkdir -p /lustre/home/ymhu/aQTL_pipeline/input/vcf/temp
	for i in {2..839}
	do
	cd /lustre/home/ymhu/aQTL_pipeline/input/vcf/temp
	#if [[ ! -f GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP_"$i" ]]
	#then
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
cd /lustre/home/ymhu/aQTL_pipeline/input
#zcat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP.vcf.gz | cut -f'$i' - > /lustre/home/ymhu/aQTL_pipeline/input/vcf/temp/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP_'$i' 
python3 /lustre/home/ymhu/aQTL_pipeline/bin/GTEx_V8_single_filter_SNPs_addcovariates.py -p /lustre/home/ymhu/aQTL_pipeline/input/vcf/temp/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_SNP_'$i' -i /lustre/home/ymhu/aQTL_pipeline/input/GTEx_covariates.v8.txt -o /lustre/home/ymhu/aQTL_pipeline/input/filter/
wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

" > /lustre/home/ymhu/aQTL_pipeline/qsub/filterVCF/filterVCF_${i}.sh
wait
qsub /lustre/home/ymhu/aQTL_pipeline/qsub/filterVCF/filterVCF_${i}.sh
echo "filterVCF_"${i}" is all finished"
#fi
wait
process=`qstat | grep "ymhu"|wc -l`
		while (( process >= 300))
		do
			echo "Current number of jobs is larger than 30"
			echo "Wait another 1 minutes"
			sleep 1m
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

run_pastefilterVCF(){
	mkdir -p /lustre/home/ymhu/aQTL_pipeline/qsub/pastefilterVCF
	mkdir -p /lustre/home/ymhu/aQTL_pipeline/output
	for i in /lustre/home/ymhu/aQTL_pipeline/Input_Covariates/Covariates_by_tissue/*.v8.covariates.txt
	do
	traitName=`echo "$i" |awk -F"/" '{print $NF;exit}'|awk -F".v8." '{print $1;exit}'`
	cd /lustre/home/ymhu/aQTL_pipeline/output
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
cd /lustre/home/ymhu/aQTL_pipeline/input/filter/
paste -d'	' id.txt ${traitName}_*.txt >/lustre/home/ymhu/aQTL_pipeline/output/${traitName}_SNPs.vcf_genotype.vcf
wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

" > /lustre/home/ymhu/aQTL_pipeline/qsub/pastefilterVCF/pastefilterVCF_${traitName}.sh

wait
qsub /lustre/home/ymhu/aQTL_pipeline/qsub/pastefilterVCF/pastefilterVCF_${traitName}.sh
echo "filterVCF_"${traitName}" is all finished"
fi
wait
process=`qstat | grep "ymhu"|wc -l`
		while (( process >= 300))
		do
			echo "Current number of jobs is larger than 30"
			echo "Wait another 1 minutes"
			sleep 1m
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
run_filterVCF2value(){
	mkdir -p /lustre/home/ymhu/aQTL_pipeline/qsub/filterVCF
	mkdir -p /lustre/home/ymhu/aQTL_pipeline/output_filter
	for i in /lustre/home/ymhu/aQTL_pipeline/output/*_SNPs.vcf_genotype.vcf
	do
	traitName=`echo "$i" |awk -F"/" '{print $NF;exit}'`
	cd /lustre/home/ymhu/aQTL_pipeline/output_filter
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
cd /lustre/home/ymhu/aQTL_pipeline/output_filter
python3 /lustre/home/ymhu/aQTL_pipeline/bin/GTEx_V8_filter_VCF_Genetype2value.py -p ${i} -o /lustre/home/ymhu/aQTL_pipeline/output_filter/${traitName}
wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

" > /lustre/home/ymhu/aQTL_pipeline/qsub/filterVCF/filterVCF_${traitName}.sh

wait
qsub /lustre/home/ymhu/aQTL_pipeline/qsub/filterVCF/filterVCF_${traitName}.sh
echo "filterVCF_"${traitName}" is all finished"
fi
wait
process=`qstat | grep "ymhu"|wc -l`
		while (( process >= 300))
		do
			echo "Current number of jobs is larger than 30"
			echo "Wait another 1 minutes"
			sleep 1m
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
run_peer(){
#Rscript bin/3UTR_QTL_peer.R ${tissue}" > ${tissue}_peer.sh

	for file in input/${tag}*combined_All_PDUIs_clean.txt
	do
	filename=${file##*/}
	tissue=${filename%_combined_All_PDUIs_clean.txt}
	cd /lustre/home/ymhu/aQTL_pipeline/peer_pdui
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
Rscript /lustre/home/ymhu/aQTL_pipeline/bin/3UTR_QTL_peer_impute.R ${tissue}
wait
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/

" > /lustre/home/ymhu/aQTL_pipeline/qsub/${tissue}_peer.sh
wait
qsub /lustre/home/ymhu/aQTL_pipeline/qsub/${tissue}_peer.sh
fi
wait
	done

#ls *expr.peer.txt|while read i; do echo -e "$i\t`awk '{print NF}' $i|head -n1`\t`awk '{print NF}' ../input/${i%.expr.peer.txt}_combined_All_PDUIs_clean.txt|head -n1`";done|less
}

run_matrixeqtl(){
#Rscript bin/run_matrixEQTL_without_plot_control_gene_exprs.R ${tissue}"> ${tissue}_v3.sh
	mkdir -p /lustre/home/ymhu/aQTL_pipeline/output_QTL
	for file in output_filter/${tag}*_SNPs.vcf_genotype.vcf
	do
	filename=${file##*/}
	tissue=${filename%_SNPs.vcf_genotype.vcf}
	cd /lustre/home/ymhu/aQTL_pipeline
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
cd /lustre/home/ymhu/aQTL_pipeline
Rscript bin/run_matrixEQTL_without_plot_covariates_final.R ${tissue}
wait
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date

wait
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
"> /lustre/home/ymhu/aQTL_pipeline/qsub/${tissue}_v3.sh

wait
qsub /lustre/home/ymhu/aQTL_pipeline/qsub/${tissue}_v3.sh
wait
		process=`qstat | grep "ymhu"|wc -l`
      	while (( process >= 300))
      	do
          	echo "Current number of jobs is larger than 300"
          	echo "Wait another 10 minutes"
          	sleep 5m
          	process=`qstat | grep "ymhu"|wc -l`
      	done
	fi
	done
}

convert_to_torus(){

#cut -f1-5 ${file} |sed 's/t.stat/t-stat/'|sed 's/p.value/p-value/' |gzip > ${file%.txt}_torus.gz"> ${tissue}_v3.sh
	for file in output_QTL/${tag}*.cis_eqtl_all_approachb.txt
	do
	filename=${file##*/}
	tissue=${filename%.cis_eqtl_all_approachb.txt}
	echo ${tissue}
	echo "#!/bin/bash
#SBATCH -A weil21_lab
#SBATCH -p standard
cut -f1-5 ${file}|sed 's/t.stat/t-stat/'|sed 's/p.value/p-value/' |gzip > ${file%.txt}_torus.gz"> ${tissue}_v3.sh
wait
sbatch ${tissue}_v3.sh

wait
		process=`squeue -u leil22|wc -l`
      	while (( process >= 100))
      	do
          	echo "Current number of jobs is larger than 30"
          	echo "Wait another 10 minutes"
          	sleep 10m
          	process=`squeue -u leil22|wc -l`
      	done
	done
}

convert_to_torus_2(){
	for file in output_QTL/${tag}*.cis_eqtl_all_approachb.txt
	do
	filename=${file##*/}
	tissue=${filename%.cis_eqtl_all_approachb.txt}
	echo ${tissue}
	echo "#!/bin/bash
#SBATCH -A weil21_lab
#SBATCH -p standard
python bin/convert_to_torus.py ${tissue}
gzip output_QTL/${tissue}.torus.trait.txt"> ${tissue}_v3.sh
wait
sbatch ${tissue}_v3.sh
wait
		process=`squeue -u leil22|wc -l`
      	while (( process >= 100))
      	do
          	echo "Current number of jobs is larger than 30"
          	echo "Wait another 10 minutes"
          	sleep 10m
          	process=`squeue -u leil22|wc -l`
      	done
	done

}


compare_with_other_covariate(){
o
	echo -e "Tissue\t3aQTL_gene_pairs\t3aQTL_gene_old_pairs\tOverlap"
	for file in output_QTL/*all.txt
	do
	filename=${file##*/}
	tissue=${filename%.cis_eqtl_all.txt}
	awk '$NF<=0.05' ${file}|cut -f1,2|cut -d"|" -f1 > tmp1.txt
	awk '$NF<=0.05' output_QTL/${tissue}.cis_eqtl_all_control_gene_exprs.txt|cut -f1,2 > tmp2.txt
	Num=`wc -l tmp1.txt|cut -d" " -f1`
	OtherNum=`wc -l tmp2.txt|cut -d" " -f1`
	Overlap=`cat tmp1.txt tmp2.txt|sort|uniq -d|wc -l`
	echo -e "${tissue}\t${Num}\t${OtherNum}\t${Overlap}"
	rm -rf -v tmp1.txt
	rm -rf -v tmp2.txt
	done

}

compare_with_other_covariate_2(){
	echo -e "Tissue\t3aQTL_gene_pairs\t3aQTL_gene_old_pairs\tOverlap"
	for file in output_QTL/*cis_eqtl_all.txt
	do
	filename=${file##*/}
	tissue=${filename%.cis_eqtl_all.txt}
	awk '$NF<=0.05' ${file}|cut -f1,2 > tmp1.txt
	line_no=`wc -l tmp1.txt|cut -d" " -f1`
	head -n ${line_no} output_QTL/${tissue}.cis_eqtl_all_approachb.txt.v2|cut -f1,2|grep -v "SNP" > tmp2.txt
	OtherNum=`wc -l tmp2.txt|cut -d" " -f1`
	Overlap=`cat tmp1.txt tmp2.txt|sort|uniq -d|wc -l`
	echo -e "${tissue}\t${line_no}\t${OtherNum}\t${Overlap}"
	rm -rf -v tmp1.txt
	rm -rf -v tmp2.txt
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



generate_qtl_plots(){
	for file in output/${tag}*_SNPs.vcf_genotype.vcf
	do
	filename=${file##*/}
	tissue=${filename%_SNPs.vcf_genotype.vcf}
	echo -e "setwd(\"/mount/weili3/Lei/Projects/2017-06-29-APA-eQTL/2018-04-12-3utr-qtl-data/output_QTL_MAF_0.1_N10/\")
	load(\"../output_QTL_MAF_0.1_N10/${tissue}.data.RDataw\")
	write.table(me\$cis\$min.pv.gene, file=paste(\"${tissue}\",\".cis.min.pv.gene.txt\",sep=\"\"))
	save(expr_subset, snps, file=paste(\"${tissue}\",\".RData\",sep=\"\"))
	save(snps, genepos, file=paste(\"${tissue}\",\".permutation.RData\",sep=\"\"))

	#cis_result <- read.table(\"${tissue}.cis_eqtl.txt\",sep=\"\\\t\",header=T)
	#cis_result <- cis_result[cis_result\$FDR<0.01,]
	#output_file_name_cis = paste0(\"${tissue}.cis_eqtl_genotype_info.txt\")
	#close( file( output_file_name_cis, open=\"w\" ) )

	#f <- function(x){
  #	snp <- as.character(x[1])
  #	gene <- as.character(x[2])
	#snp_info = snps\$FindRow(snp)\$row
  #	s1 = as.numeric(snp_info)
	#write.table(snp_info, file=output_file_name_cis, append=TRUE, sep=\"\\\t\")
	#}

	#pdf(\"${tissue}.cis.eQTL.pdf\")
	#apply(cis_result, 1, f)
	#dev.off()
	" > ${tissue}.R

	echo "#!/bin/bash
	#PBS -m e
	#PBS -M lei.li@bcm.edu
	#PBS -k oe
	#PBS -l nodes=1:ppn=8
	#PBS -l vmem=20gb
	#PBS -l walltime=40:00:00
	module load R/3.4.0
	cd \$PBS_O_WORKDIR
	Rscript ${tissue}.R" > ${tissue}_plot.sh

	done

}

generate_trans_qtl_plots(){
	for file in output/${tag}*_SNPs.vcf_genotype.vcf
	do
	filename=${file##*/}
	tissue=${filename%_SNPs.vcf_genotype.vcf}
	echo ${tissue}
	echo -e "setwd(\"/mount/weili3/Lei/Projects/2017-06-29-APA-eQTL/2018-04-12-3utr-qtl-data/output_QTL\")
	load(\"${tissue}.data.RDataw\")
	trans_result <- read.table(\"${tissue}.trans_eqtl.txt\",sep=\"\\\t\",header=T)
	trans_result <- trans_result[trans_result\$FDR<1e-4,]
	output_file_name_trans = paste0(\"${tissue}.trans_eqtl_genotype_info.txt\")
	close( file( output_file_name_trans, open=\"w\" ) )

	f <- function(x){
  	snp <- as.character(x[1])
  	gene <- as.character(x[2])
  	e1 = as.numeric(expr_subset[which(rownames(expr_subset)==gene),])
	  snp_info = snps\$FindRow(snp)\$row
  	s1 = as.numeric(snp_info)
	  write.table(snp_info, file=output_file_name_trans, append=TRUE, sep=\"\\\t\")
  	lm1 = lm(e1 ~ s1)
  	plot(e1 ~ jitter(s1),col=(s1+1),xaxt=\"n\",xlab=\"Genotype\",ylab=\"PDUI\",main=paste(snp,gene,sep=\"|||\"))
  	axis(1,at=c(0:2),labels=c(\"REF\",\"HET\",\"ALT\"))
	}

	pdf(\"${tissue}.trans.eQTL.pdf\")
	apply(trans_result, 1, f)
	dev.off()
	" > ${tissue}.R

	echo "#!/bin/bash
	#PBS -m e
	#PBS -M lei.li@bcm.edu
	#PBS -k oe
	#PBS -l nodes=1:ppn=8
	#PBS -l vmem=20gb
	#PBS -l walltime=30:00:00
	module load R/3.4.0
	cd \$PBS_O_WORKDIR
	Rscript ${tissue}.R" > ${tissue}_plot.sh

	done
}

get_sig_qtl_list(){

	for file in output_QTL/*cis_eqtl.txt
	do
	awk '$NF<0.05' ${file}|cut -f1|sort|uniq|grep -v "SNP" > ${file%.txt}_sig_id.txt
	done

}

generate_subset_genotype(){
	for file in output/*_SNPs.vcf_genotype.vcf
	do
	filename=${file##*/}
	tissue=${filename%_SNPs.vcf_genotype.vcf}
	echo "#!/bin/bash
#PBS -m e
#PBS -M lei.li@bcm.edu
#PBS -k oe
#PBS -l nodes=1:ppn=8
#PBS -l vmem=80gb
#PBS -l walltime=10:00:00
module load R/3.4.3
cd \$PBS_O_WORKDIR
Rscript bin/extract_genotype.r ${tissue}" > ${tissue}_genotype.sh
	done

}

run_permutation_test(){
	for Rdata in output_QTL/*.permutation.RData
	do
	filename=${Rdata##*/}
	tissue=${filename%.permutation.RData}
	echo ${tissue}
			for i in `seq 111 114`
			do
				echo "#!/bin/bash
#$ -q gpu*,free*,pub*
#$ -pe openmp 8-64
module load R/3.4.1
Rscript bin/_eQTL_permutation_minP.R ${i} ${Rdata} peer_pdui/${tissue}.expr.peer.txt ${tissue} " > ${tissue}_${i}_genotype.sh
			qsub -ckpt restart ${tissue}_${i}_genotype.sh
				process=`qstat|grep "leil22"|wc -l`
      	while (( process >= 250))
      	do
          	echo "Current number of jobs is larger than 30"
          	echo "Wait another 30 minutes"
          	sleep 10m
          	process=`qstat|grep "leil22"|wc -l`
      	done
			done
		done
}

run_vcf_filtering(){
	for file in /mount/weili3/Lei/Projects/2017-06-29-APA-eQTL/2018-04-12-3utr-qtl-data/output/*_SNPs.vcf.gz
	do
	filename=${file##*/}_vcf
	tbi_file=output/${filename%_SNPs.vcf.gz_vcf}_SNPs.filtered.vcf.gz.tbi
	if [ ! -f ${tbi_file} ]
	then
	echo $filename
	echo "#!/bin/bash" > ${filename}.sh
	echo "#PBS -m e" >> ${filename}.sh
	echo "#PBS -M lei.li.bioinfo@gmail.com" >> ${filename}.sh
	echo "#PBS -k oe" >> ${filename}.sh
	echo "#PBS -l nodes=1:ppn=1" >> ${filename}.sh
	echo "#PBS -l vmem=20gb" >> ${filename}.sh
	echo "#PBS -l walltime=20:00:00" >> ${filename}.sh
	#echo "module load bcftools/1.3" >> ${filename}.sh
	#echo "cd \$PBS_O_WORKDIR" >> ${filename}.sh
	echo "bcftools view -q 0.01:minor -o ${file%.vcf.gz}.filtered.vcf.gz -O z ${file}" >> ${filename}.sh
	echo "tabix -p vcf ${file%.vcf.gz}.filtered.vcf.gz" >> ${filename}.sh
	#echo "tabix -p vcf ${file}" >> ${filename}.sh
	fi
	done

}

get_final_qtl(){
	for Rdata in output_QTL/*cis.min.pv.gene.txt
	do
		filename=${Rdata##*/}
		tissue=${filename%.cis.min.pv.gene.txt}
		echo $tissue
		echo "#!/bin/bash
#$ -q gpu*,free*,pub*
#$ -pe openmp 8-64
module load R/3.4.1
Rscript bin/permutation_analysis.R ${tissue}" > ${tissue}_v.sh
		qsub -ckpt restart ${tissue}_v.sh
		process=`qstat|grep "leil22"|wc -l`
      	while (( process >= 100))
      	do
          	echo "Current number of jobs is larger than 30"
          	echo "Wait another 30 minutes"
          	sleep 1m
          	process=`qstat|grep "leil22"|wc -l`
      	done
	done

}

get_final_qtl_v2(){
	# This is based on already generated APAgenes by last model "get_final_qtl"
	for Rdata in output_QTL/*cis_eqtl_all.txt
	do
		filename=${Rdata##*/}
		tissue=${filename%.cis_eqtl_all.txt}
		echo ${tissue}
		awk '$NF<0.05' ${Rdata} > ${Rdata%_all.txt}_FDR0.05.txt
		agene="../2018-04-12-3utr-qtl-data/output_QTL_MAF_0.1_N10/${tissue}.cis.agene.txt"
		cat ${agene}|while read j; do grep $j ${Rdata%_all.txt}_FDR0.05.txt; done > output_QTL/${tissue}.cis.3QTL.FDR.05.xls
		rm ${Rdata%_all.txt}_FDR0.05.txt
	done
}

get_snp_count(){
	for Rdata in output/${tag}*_SNPs.vcf_genotype.vcf
	do
	filename=${Rdata##*/}
        tissue=${filename%_SNPs.vcf_genotype.vcf}
	echo ${tissue}
	if [ ! -f output_QTL_MAF_0.1_N10/${tissue}_all_snp_count.txt ]; then
	cut -f2 output_QTL_MAF_0.1_N10/${tissue}.cis_eqtl_all.txt |sort|uniq -c > output_QTL_MAF_0.1_N10/${tissue}_all_snp_count.txt
	cut -f2 output_QTL_MAF_0.1_N10/${tissue}.cis.3QTL.FDR.05.xls|sort|uniq -c > output_QTL_MAF_0.1_N10/${tissue}_signif_count.txt
	fi
	done


}

eGene_filter(){

	for Rdata in output_QTL_MAF_0.1_N10/*cis.min.pv.gene.txt
	do
	filename=${Rdata##*/}
	tissue=${filename%.cis.min.pv.gene.txt}
	echo ${tissue}
	python bin/GTEx_eGene_filter.py ${tissue}
	done
}


add_rsid(){
	python bin/find_rs.py /mount/weili2/Lei/SNPs_RSnames/ sig_3QTL/

}

plot_specific_genes(){
	#tissue="Muscle_Skeletal"
	#Rscript bin/QTL_plot.R -t ${tissue} -s "chr19_8501092_C_T" -g "NM_001005415|MARCH2|chr19|+"
	#Rscript bin/QTL_plot.R -t ${tissue} -s "chr5_175370893_C_T" -g "NM_032361|THOC3|chr5|-"

	#chr22_21936152_C_G      NR_046082|UBE2L3|chr22|+ Whole_Blood
	#Rscript bin/QTL_plot.R -t "Whole_Blood" -s "chr22_21936152_C_G" -g "NR_046082|UBE2L3|chr22|+"


	#Rscript bin/QTL_plot.R -t "Adipose_Subcutaneous" -s "chr6_31238147_C_G" -g "NM_002117|HLA-C|chr6|-"

	#Rscript bin/QTL_plot.R -t "Adipose_Subcutaneous" -s "chr6_29913639_A_G" -g "NM_001242758|HLA-A|chr6|+"

	for Rdata in output_QTL/*cis.min.pv.gene.txt
        do
        filename=${Rdata##*/}
        tissue=${filename%.cis.min.pv.gene.txt}
        echo ${tissue}
        #Rscript bin/QTL_plot.R -t ${tissue} -s "chr7_128589427_G_A" -g "NM_001098630|IRF5|chr7|+"
        #Rscript bin/QTL_plot.R -t ${tissue} -s "chr22_21936152_C_G" -g "NR_046082|UBE2L3|chr22|+"
				Rscript bin/QTL_plot.R -t ${tissue} -s "chr6_29913639_A_G" -g "NM_001242758|HLA-A|chr6|+"
	done

}

find_lead_3aQTL(){
		for qtlf in output_QTL/*cis.3QTL.FDR.05.xls
		do
			echo ${qtlf}
			filename=${qtlf##*/}
			tissue=${filename%.cis.3QTL.FDR.05.xls}
			sort -k2,2 -k6,6n ${qtlf} > ${qtlf%.xls}.sort.xls
			~/anaconda2/bin/python bin/find_best_v2.py ${qtlf%.xls}.sort.xls output_QTL/${tissue}.best.3QTL.txt
		done

}

find_lead_3aQTL_2(){
		#for qtlf in output_QTL/*.cis_eqtl_all.txt
		#for qtlf in output_QTL/*.cis_eqtl_all_approachb.txt
		for qtlf in output_QTL/*.torus.trait.txt.gz
		do
			echo ${qtlf}
			filename=${qtlf##*/}
			#tissue=${filename%.cis_eqtl_all.txt}
			#tissue=${filename%.cis_eqtl_all_approachb.txt}
                        tissue=${filename%.torus.trait.txt.gz}
echo "#!/bin/bash
#SBATCH -A weil21_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#awk '\$NF<0.05' ${qtlf}| sort -k2,2 -k6,6g  > ${qtlf%.txt}.sort.xls
gunzip ${qtlf}
awk '\$NF<1e-4' ${qtlf%.gz}|sort -k2,2 -k5,5g > ${qtlf%.txt.gz}.sort.xls 
python bin/find_best_v2.py ${qtlf%.txt.gz}.sort.xls output_QTL/${tissue}.best.LDSR_non1K.3aQTL.txt
gzip ${qtlf%.gz}" > ${tissue}_lead.sh
#python bin/find_best_v2.py ${qtlf%.txt}.sort.xls output_QTL/${tissue}.best.isoQTL.txt" > ${tissue}_lead.sh
wait
sbatch ${tissue}_lead.sh
wait
		done
}

compare_best_variants(){
	for qtlf in output_QTL/*.cis_eqtl_all_approachb.txt
		do
			echo ${qtlf}
			filename=${qtlf##*/}
			#tissue=${filename%.cis_eqtl_all.txt}
			tissue=${filename%.cis_eqtl_all_approachb.txt}
			python bin/Compare_isoQTL_aQTL.py ${tissue} > ${tissue}_compare.txt
		done
}

3QTLs_to_publish(){
	mkdir -p files_to_publish
	#for file in output_QTL/*3QTL.FDR.05.xls
	#do
	#	filename=${file##*/}
	#	cut -f1,2,3,4,5 ${file} > files_to_publish/${filename%3QTL.FDR.05.xls}3aQTL.txt
	#done

	for file in output_QTL/Pituitary*.cis_eqtl_all_control_gene_exprs.txt
	do
		filename=${file##*/}
		tissue=${filename%.cis_eqtl_all_control_gene_exprs.txt}
		echo -e ""SNP"\t"gene"\t"beta"\t"t-stat"\t"p-value"" > files_to_publish/${tissue}_3aQTL_expr_controlled.txt
		awk '$NF<0.05' $file|cut -f1,2,3,4,5 >> files_to_publish/${tissue}_3aQTL_expr_controlled.txt
	done


}

convert_3QTL_to_bed(){
	mkdir -p bed_files
	for file in output/${tag}*_SNPs.vcf_genotype.vcf
	do
	filename=${file##*/}
	tissue=${filename%_SNPs.vcf_genotype.vcf}
	python bin/Convert_to_bed_format.py ${tissue}
	done

}

compare_qtl_differenes(){
		#echo -e ""Tissue"\t"Covariant_adjusted"\t"Old"\t"Overlap_Number"" > tissue_pairs_number.txt
		#for file in output_QTL/*3QTL.FDR.05.xls
		#do
		#filename=${file##*/}
		#tissue=${filename%.cis.3QTL.FDR.05.xls}
		#prev_dat="../2018-04-12-3utr-qtl-data/output_QTL_MAF_0.1_N10/${tissue}.cis.3QTL.FDR.05.xls"

		#cut -f1,2 ${file}|grep -v "gene" > tmp1.txt
		#cut -f1,2 ${prev_dat}|grep -v "gene" > tmp2.txt
		#overlap_nu=`cat tmp1.txt tmp2.txt|sort|uniq -d|wc -l`
		#echo -e "${tissue}\t`wc -l tmp1.txt|cut -d " " -f1`\t`wc -l tmp2.txt|cut -d " " -f1`\t${overlap_nu}" >> tissue_pairs_number.txt
		#rm tmp1.txt tmp2.txt
		#done


		#echo -e ""Tissue"\t"Covariant_adjusted"\t"Old"\t"Overlap_Number"" > tissue_pairs_number_QTL_all.txt
                #for file in output_QTL/*.cis_eqtl_all.txt
                #do
                #filename=${file##*/}
                #tissue=${filename%.cis_eqtl_all.txt}
                #prev_dat="../2018-04-12-3utr-qtl-data/output_QTL_MAF_0.1_N10/${tissue}.cis_eqtl_all.txt"

                #cut -f1,2 ${file} > tmp1.txt
                #cut -f1,2 ${prev_dat}> tmp2.txt
                #overlap_nu=`cat tmp1.txt tmp2.txt|sort|uniq -d|wc -l`
                #echo -e "${tissue}\t`wc -l tmp1.txt|cut -d " " -f1`\t`wc -l tmp2.txt|cut -d " " -f1`\t${overlap_nu}" >> tissue_pairs_number_QTL_all.txt
                #rm tmp1.txt tmp2.txt
                #done


		#echo -e ""Tissue"\t"Covariant_adjusted"\t"Old"\t"Overlap_Number"" > tissue_pairs_number_QTL_test.txt
		#for file in output_QTL_sex_pc/*.cis_eqtl_all.txt
		#do
		#filename=${file##*/}
		#tissue=${filename%.cis_eqtl_all.txt}
		#prev_dat="../2018-04-12-3utr-qtl-data/output_QTL_MAF_0.1_N10/${tissue}.cis_eqtl_all.txt"
		#prev_dat="output_QTL_sex/${tissue}.cis_eqtl_all.txt"
		#awk '$NF<0.05' ${file}|cut -f1,2 > tmp1.txt
		#awk '$NF<0.05' ${prev_dat}|cut -f1,2 > tmp2.txt
		#overlap_nu=`cat tmp1.txt tmp2.txt|sort|uniq -d|wc -l`
		#echo -e "${tissue}\t`wc -l tmp1.txt|cut -d " " -f1`\t`wc -l tmp2.txt|cut -d " " -f1`\t${overlap_nu}" >> tissue_pairs_number_QTL_test.txt
		#rm tmp1.txt tmp2.txt
		#done


		for qtlf in output_QTL/*cis_eqtl*.txt
		do
			echo ${qtlf}
			filename=${qtlf##*/}
			awk '$NF<0.05' ${qtlf} > ${qtlf%.txt}.tmp
			sort -k2,2 -k6,6n ${qtlf%.txt}.tmp > ${qtlf%.txt}.sort.tmp
			python bin/find_best_v2.py ${qtlf%.txt}.sort.tmp ${qtlf%.txt}.sort.best.tmp
			rm ${qtlf%.txt}.tmp
			rm ${qtlf%.txt}.sort.tmp
		done

		#for qtlf in output_QTL/*cis_eqtl_all.txt
		#do
		#	echo ${qtlf}
		#	filename=${qtlf##*/}
		#	tissue=${filename%.cis_eqtl_all.txt}
		#	output=output_QTL/${tissue}.cis_eqtl_qtl_pairs.txt
		#	if [ ! -f ${output} ]
		#	then
		#	awk '$NF<0.05' ${qtlf} > ${qtlf%all.txt}fdr0.05.txt
		#	cut -f1,2 ${qtlf%all.txt}fdr0.05.txt > ${qtlf%all.txt}qtl_pairs.txt
		#	fi

			#cut -f2 output_QTL_sex/${tissue}.cis.3QTL.FDR.05.xls |sort|uniq|while read gene;do grep $gene ${qtlf%all.txt}fdr0.05.txt ;done > ${qtlf%.cis_eqtl_all.txt}.cis.3QTL.FDR.05.xls
			#~/anaconda2/bin/python bin/Compare_qtl_pairs.py ${tissue}

		#	~/anaconda2/bin/python bin/find_best_v2.py ${qtlf%.cis_eqtl_all.txt}.cis.3QTL.FDR.05.xls ${qtlf%.cis_eqtl_all.txt}.best.3QTL.txt
		#	awk '$NF<1e-06' ${qtlf%.cis_eqtl_all.txt}.best.3QTL.txt|cut -f1,2 |sed 's/\s/\t/g' > ${qtlf%.cis_eqtl_all.txt}.best.3QTL.tmp

		#		cut -f2 output_QTL_sex/${tissue}.cis.3QTL.FDR.05.xls |sort|uniq|while read gene; do grep $gene output_QTL_sex/${tissue}.cis.3QTL.FDR.05.xls; done > output_QTL_sex/${tissue}.cis.3QTL.FDR.05.xls.tmp
		#	~/anaconda2/bin/python bin/find_best_v2.py output_QTL_sex/${tissue}.cis.3QTL.FDR.05.xls.tmp output_QTL_sex/${tissue}.best.3QTL.txt

		#	awk '$NF<1e-06' output_QTL_sex/${tissue}.best.3QTL.txt|cut -f1,2 |sed 's/\s/\t/g' > output_QTL_sex/${tissue}.best.3QTL.tmp
		#	~/anaconda2/bin/python bin/Compare_qtl_pairs.py ${tissue}
		#done
}

main

