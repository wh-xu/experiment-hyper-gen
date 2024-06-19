import sys, os, glob

def get_largest_fna_files(folder_path, num_files):
    files = glob.glob(os.path.join(folder_path, "*.fna"))
    file_sizes = {file: os.path.getsize(file) for file in files}
    files_sorted = sorted(file_sizes.items(), key=lambda x: x[1], reverse=True)
    
    files = [os.path.basename(i[0]) for i in files_sorted[:num_files]]
    for file in files:
        print(file)
        

get_largest_fna_files(sys.argv[1], int(sys.argv[2]))
