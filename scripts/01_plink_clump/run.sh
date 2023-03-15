main(){
    global_parameters
    run_plink_clump
    run_finemapping
}

global_parameters(){
    CUR_DIR="/lustre/home/hchen/cancer-GWAS"

}

run_plink_clump(){
  cd $CUR_DIR/01plink_clump/output
  module load anaconda/anaconda3-2019.10-py37
  python $CUR_DIR/01plink_clump/bin/run_plink_clump.py -p $CUR_DIR/01plink_clump/bin/allGWAS_number_files.txt >run_plink_clump.log 2>run_plink_clump.err
}

run_finemapping(){
  cd $CUR_DIR/output/CAUSALdb
  module load anaconda/anaconda3-2019.10-py37
  for file in $CUR_DIR/output/CAUSALdb/input/*.txt
  do
      name=`echo "$file" | awk -F "/" '{print $NF }'|awk -F"." '{print $1;exit}'`
        cd $CUR_DIR/output/CAUSALdb/input
        if [[ ! -f GLGC_"$name"_result.txt ]]
        then
        cp $file $CUR_DIR/output/CAUSALdb/input
        mv "$name".txt $CUR_DIR/output/CAUSALdb/input/GLGC_${name}_result.txt  ### unity the name of GLGC_CG0XXX_result.txt
        fi
done
  cd $CUR_DIR/output/CAUSALdb
  python $CUR_DIR/output/CAUSALdb/run_fine_mapping.py -p $CUR_DIR/bin/FinnGen_r5_numberfile.txt >run_fine_mapping.log 2>run_fine_mapping.err &     
}


main
