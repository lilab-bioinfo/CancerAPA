main(){ 
  compute_weights
  generate_summary
  run_3aTWAS_analysis
  conditional_analysis
}

compute_weights(){
for tiss in `ls tissue.list`
  do
    WORK_DIR=`pwd`
    INPUT_DIR=$WORK_DIR/input
    OUTPUT_DIR=$WORK_DIR/output
    PRE="${tiss}"
    
    mkdir -p $WORK_DIR/src/slurm/${tiss}
    
    PRE_GECP="${PRE}_combined_All_PDUIs_clean.txt"
    BATCH_START=1
    BATCH_END=`less ${INPUT_DIR}/APA_matrix/$PRE_GERP}|wc -l`
    NR="${BATCH_START}_${BATCH_END}"
    
    
  
  
done
}

generate_summary(){





}

run_3aTWAS_analysis(){





}

contidional_analysis(){




}


main
