OVERWRITE=true
TOPK=100

for ds in D1 D2 D3 D5 GTDB_r207;
do
    DS_PATH=../dna-dataset/$ds
    NEW_PATH=${DS_PATH}_$TOPK
    echo $NEW_PATH

    if $OVERWRITE; then
        mkdir -p $NEW_PATH
        python get_largest_fna_files.py $DS_PATH $TOPK >$NEW_PATH/fna_files.txt

        for file in $(cat $NEW_PATH/fna_files.txt); do
            echo "Moving $file to $NEW_PATH"
            cp $DS_PATH/$file $NEW_PATH
        done
    fi

    realpath $NEW_PATH/*.fna > $NEW_PATH/ref_files.txt
    cp $NEW_PATH/ref_files.txt $NEW_PATH/query_files.txt
done
