main(){
        #run_get_bedgrach
        #create_subtissues
        #check_wig
        #dapars2_analysis
        run_peer
}

run_get_bedgraph(){
        for file in Data/*.bigWig
        do
        filename=${file##*/}
        tissue=${filename%.Aligned.sortedByCoord.out.patched.md.bigWig}
        echo ${tissue}
        if [ ! -f ${file%.out.patched.md.bigWig}.wig ]
        then
        echo "#!/bin/bash
#SBATCH -A weil21_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
bigWigToBedGraph $file ${file%.out.patched.md.bigWig}.wig
bzip2 $file" > scripts/${tissue}_peer.sh
wait
sbatch scripts/${tissue}_peer.sh
wait
        process=`squeue -u leil22|wc -l`
        while (( process >= 300))
        do
                echo "Current number of jobs is larger than 30"
                echo "Wait another 10 minutes"
                sleep 1m
                process=`squeue -u leil22|wc -l`
        done
    fi
        done

#ls *expr.peer.txt|while read i; do echo -e "$i\t`awk '{print NF}' $i|head -n1`\t`awk '{print NF}' ../input/${i%.expr.peer.txt}_combined_All_PDUIs_clean.txt|head -n1`";done|less
}

check_wig(){
        for file in Data/*.wig
        do
        filename=${file##*/}
        tissue=${filename%.wig}
        echo ${tissue}
        echo "#!/bin/bash
#SBATCH -A weil21_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#cut -f1 ${file}|uniq > ${file%.wig}.txt
wc -l ${file%.wig}.txt > ${file%.wig}.count" > scripts/${tissue}_peer.sh
wait
sbatch scripts/${tissue}_peer.sh
wait
        process=`squeue -u leil22|wc -l`
        while (( process >= 600))
        do
                echo "Current number of jobs is larger than 30"
                echo "Wait another 10 minutes"
                sleep 1m
                process=`squeue -u leil22|wc -l`
        done
        done


}

create_subtissues(){
        for file in ~/data/GTEx/subTissues/*.sampleID
        do
                filename=${file##*/}
                tissue=${filename%.sampleID}
                mkdir -p ${tissue}/
                cp ${file} ${tissue}/${tissue}
        done


}

dapars_analysis(){
files=(chrX chrY chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 
chr19 chr20 chr21 chr22)
for (( i = 0 ; i < ${#files[@]} ; i++ ))
do
        filename=${files[$i]}
        echo "#!/bin/bash
#SBATCH -A weil21_lab
#SBATCH -p highmem
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=3
module load R/3.6.2
module load bedtools2/2.29.2
~/anaconda2/bin/python ~/GTEx/src/DaPars2_Multi_Sample_Multi_Chr.py mapping_joint_Analysis_DaPars_each_chr_configure.txt ${files[$i]}" > ${filename}_dapars.sh
wait
sbatch ${filename}_dapars.sh
wait
done



}

run_peer(){
        #sort -k1,1 -u
        for file in PDUIs/*_clean.txt
        do
        filename=${file##*/}
        tissue=${filename%_combined_All_PDUIs_clean.txt}
        if [ ! -f peer_pdui/${tissue}.expr.peer.txt ]
        then
        echo ${tissue}
        echo "#!/bin/bash
#SBATCH -A weil21_lab
#SBATCH -p standard
#SBATCH --nodes=4
#SBATCH --cpus-per-task=4
module load R/4.0.2
Rscript src/3UTR_QTL_peer.R ${tissue}" > ${tissue}_peer.sh
wait
sbatch ${tissue}_peer.sh
wait
        fi
        done

#ls *expr.peer.txt|while read i; do echo -e "$i\t`awk '{print NF}' $i|head -n1`\t`awk '{print NF}' ../input/${i%.expr.peer.txt}_combined_All_PDUIs_clean.txt|head -n1`";done|less       

}

main
