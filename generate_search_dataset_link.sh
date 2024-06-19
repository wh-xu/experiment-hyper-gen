OVERWRITE=true

for ds in D1 D2 D3 D5;
do
    DS_PATH=../dna-dataset/$ds
    NEW_PATH=../dna-dataset/${ds}_all
    echo "Generating ref/query list for $NEW_PATH"
    # ln -s $DS_PATH $NEW_PATH
    # realpath $NEW_PATH/*.fna > $NEW_PATH/ref_files.txt

    if [ $ds == D1 ]
    then
        realpath ../dna-dataset/D3/Escherichia_coli_str__K_12_substr__W3110_NC_007779.LargeContigs.fna > $NEW_PATH/query_files.txt
    fi

    if [ $ds == D2 ]
    then
        realpath $DS_PATH/Bacillus_anthracis_52_G_NZ_CM002395.LargeContigs.fna > $NEW_PATH/query_files.txt
    fi

    if [ $ds == D3 ]
    then
        realpath $DS_PATH/Escherichia_coli_0_1288_GCA_000303255.LargeContigs.fna > $NEW_PATH/query_files.txt
    fi

    if [ $ds == D5 ]
    then
        realpath ../dna-dataset/D3/Escherichia_coli_str__K_12_substr__W3110_NC_007779.LargeContigs.fna > $NEW_PATH/query_files.txt
    fi

done
