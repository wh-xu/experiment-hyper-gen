target_dir=$1
echo $target_dir
for d in $target_dir/*/ ; 
do
    # Loop through each subdirectory
    # and count the number of words in each file
#    (cd "$d" && wc -w *.txt);
    (cd "$d" && echo $(pwd));
done
