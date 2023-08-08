main(){
    global_parameters
    run_plink_clump
}

global_parameters(){
    CUR_DIR=`pwd`

}

run_plink_clump(){
  cd $CUR_DIR/01plink_clump/output
  module load anaconda/anaconda3-2019.10-py37
  python $CUR_DIR/01plink_clump/bin/run_plink_clump.py -p $CUR_DIR/01plink_clump/bin/allGWAS_number_files.txt >run_plink_clump.log 2>run_plink_clump.err
}

main
