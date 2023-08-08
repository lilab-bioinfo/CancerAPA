main(){
    global_parameters
    run_finemapping
}

global_parameters(){
    CUR_DIR=`pwd`

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
        mv "$name".txt $CUR_DIR/CAUSALdb/input/GLGC_${name}_result.txt  ### unity the name of GLGC_CG0XXX_result.txt
        fi
done
  cd $CUR_DIR/output/CAUSALdb
  python $CUR_DIR/CAUSALdb/run_fine_mapping.py -p $CUR_DIR/CAUSALdb/fine_mapping_input.txt >run_fine_mapping.log 2>run_fine_mapping.err &     
}


main
