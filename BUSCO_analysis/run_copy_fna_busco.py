import os
import glob
import tqdm

file = "BUSCO_fna_list.txt"

with open(file) as fp:
    fna_file = fp.read().splitlines()


path_dataset = [
    "/project/weihong/dna-dataset/D1_all",
    "/project/weihong/dna-dataset/D2_all",
    "/project/weihong/dna-dataset/D3_all",
    "/project/weihong/dna-dataset/D5_all",
    "/project/weihong/dna-dataset/GTDB_r207_all",
]

all_fna_files = []
for p in path_dataset:
    all_fna_files.extend(glob.glob(os.path.join(p, "*.fna")))

print(len(all_fna_files))

os.makedirs("fna_files", exist_ok=True)

for i, f in tqdm.tqdm(enumerate(fna_file)):
    for r in all_fna_files:
        if f in r:
            os.system(f"cp {r} fna_files/")

# busco -c 16 -m genome -i ./fna_files --out_path ./results
