K="21"
THREADS=16
SKETCH_SIZE="512 1024 2048 4096 8192"



for tool in hd
do
    METHOD="t1ha2"
    HD_D="1024 2048 4096 8192"
    HD_S="1000 1500 2000"
    DEVICE="gpu"

    for ds in D1 D2 D3 D5 GTDB_r207
    do
        SUFFIX=${ds}_all
        REF_LIST=../dna-dataset/$SUFFIX/ref_files.txt
        QUERY_LIST=../dna-dataset/$SUFFIX/query_files.txt

        echo "${SUFFIX}"
        python main_ani_experiment.py --exp_type search -p $THREADS -k $K -s $SKETCH_SIZE --hd_scaled $HD_S --hd_dim $HD_D --device $DEVICE --method $METHOD -query_list  $QUERY_LIST -ref_list $REF_LIST --only_type $tool -suffix $SUFFIX > ./exp_log/exp_db_${SUFFIX}_${tool}.log
    done
done



for tool in skani fastani bindash dashing2 mash 
do
    METHOD="bagminhash"

    for ds in D1 D2 D3 D5 GTDB_r207
    do
        SUFFIX=${ds}_all
        REF_LIST=../dna-dataset/$SUFFIX/ref_files.txt
        QUERY_LIST=../dna-dataset/$SUFFIX/query_files.txt

        echo "${SUFFIX}"
        python main_ani_experiment.py --exp_type search -p $THREADS -k $K -s $SKETCH_SIZE --hd_scaled $HD_S --hd_dim $HD_D --device $DEVICE --method $METHOD -query_list  $QUERY_LIST -ref_list $REF_LIST --only_type $tool -suffix $SUFFIX > ./exp_log/exp_db_${SUFFIX}_${tool}.log
    done
done

