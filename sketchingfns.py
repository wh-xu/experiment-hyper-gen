import sys
import os
import pandas as pd
import numpy as np
import subprocess
from joblib import Parallel, delayed
import glob
from timeit import default_timer as time
import random
import tqdm
import glob


def getrstr():
    return "".join(map(str, (random.choice(range(10)) for x in range(20))))


def check_call(*args, **kwargs):
    trynum = 0
    while 1:
        try:
            return subprocess.check_output(*args, **kwargs)
        except subprocess.CalledProcessError as e:
            trynum += 1
            print("Failed call once; trying for %dth time" % trynum, file=sys.stderr)
            if trynum == 3:
                raise
            else:
                continue


def repeat_x(func, ntimes, *args, **kwargs):
    results = [func(*args, **kwargs) for i in range(ntimes)]
    return np.median(results)


def func_run_anim(query_files, ref_files, path_out, threads):
    # Parse NUCmer delta file to get total alignment length and total sim_errors
    def parse_delta(filename):
        """Return (alignment length, similarity errors) tuple from passed .delta.

        :param filename:  Path, path to the input .delta file

        Extracts the aligned length and number of similarity errors for each
        aligned uniquely-matched region, and returns the cumulative total for
        each as a tuple.

        Similarity errors are defined in the .delta file spec (see below) as
        non-positive match scores. For NUCmer output, this is identical to the
        number of errors (non-identities and indels).

        Delta file format has seven numbers in the lines of interest:
        see http://mummer.sourceforge.net/manual/ for specification

        - start on query
        - end on query
        - start on target
        - end on target
        - error count (non-identical, plus indels)
        - similarity errors (non-positive match scores)
            [NOTE: with PROmer this is equal to error count]
        - stop codons (always zero for nucmer)

        To calculate alignment length, we take the length of the aligned region of
        the reference (no gaps), and process the delta information. This takes the
        form of one value per line, following the header sequence. Positive values
        indicate an insertion in the reference; negative values a deletion in the
        reference (i.e. an insertion in the query). The total length of the alignment
        is then:

        reference_length + insertions - deletions

        For example:

        A = ABCDACBDCAC$
        B = BCCDACDCAC$
        Delta = (1, -3, 4, 0)
        A = ABC.DACBDCAC$
        B = .BCCDAC.DCAC$

        A is the reference and has length 11. There are two insertions (positive delta),
        and one deletion (negative delta). Alignment length is then 11 + 1 = 12.
        """
        in_aln, aln_length, sim_errors = False, 0, 0
        for line in [_.strip().split() for _ in open(filename, "r").readlines()]:
            if line[0] == "NUCMER" or line[0].startswith(">"):  # Skip headers
                continue
            # Lines with seven columns are alignment region headers:
            if len(line) == 7:
                aln_length += abs(int(line[1]) - int(line[0])) + 1  # reference length
                sim_errors += int(line[4])  # count of non-identities and indels
                in_aln = True
            # Lines with a single column (following a header) report numbers of symbols
            # until next insertion (+ve) or deletion (-ve) in the reference; one line per
            # insertion/deletion; the alignment always ends with 0
            if in_aln and line[0].startswith("0"):
                in_aln = False
            elif in_aln:
                # Add one to the alignment length for each reference insertion; subtract
                # one for each deletion
                val = int(line[0])
                if val < 1:  # deletion in reference
                    aln_length += 1
                elif val == 0:  # ends the alignment entry
                    in_aln = False
        return aln_length, sim_errors

    def generate_anim_job_list(query_files, ref_files, path_out):
        job_list = []
        for i, query in enumerate(query_files):
            for j, ref in enumerate(ref_files):
                prefix = (
                    os.path.splitext(os.path.basename(query))[0]
                    + "_vs_"
                    + os.path.splitext(os.path.basename(ref))[0]
                )
                pathf = os.path.join(path_out, "nucmer_output", prefix)
                cmd = f"nucmer --mum -p {pathf} {query} {ref}"
                job_list.append(cmd)

        return job_list

    def execute_cmd(cmd):
        return subprocess.check_call(
            cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )

    os.makedirs(os.path.join(path_out, "nucmer_output"), exist_ok=True)

    startt = time()
    job_list = generate_anim_job_list(query_files, ref_files, path_out)
    Parallel(n_jobs=threads)(
        delayed(execute_cmd)(i) for i in tqdm.tqdm(job_list, leave=False)
    )
    run_time = time() - startt

    deltafiles = glob.glob(os.path.join(path_out, "nucmer_output", "*.delta"))
    df = pd.DataFrame(columns=["query", "reference", "ani"])

    for deltafile in tqdm.tqdm(deltafiles, leave=False):
        qname, rname = deltafile[:-6].split("_vs_")
        tot_length, tot_sim_error = parse_delta(deltafile)
        perc_id = 1 - float(tot_sim_error) / tot_length if tot_length != 0 else 0
        df.loc[len(df)] = qname, rname, perc_id
    df.to_csv(os.path.join(path_out, "pw_ani.csv"), index=False)

    return run_time


def func_run_fastani(
    query_files, ref_files, path_out, file_out, threads, exe_path="./tools/fastANI"
):
    startt = time()
    check_call(
        f"{exe_path} -t {threads} --ql {query_files} --rl {ref_files} -o {os.path.join(path_out, file_out)}",
        shell=True,
    )
    return time() - startt


def func_mash_sketch(file_out, file_list, k, size, threads, exe_path="./tools/mash"):
    starttm = time()
    check_call(
        f"{exe_path} sketch -s {size} -k {k} -o {file_out} -l {file_list} -p {threads} 2> /tmp/tmp.dat",
        shell=True,
    )
    return time() - starttm


def func_mash_dist(sketch_file, threads, dist_file, exe_path="./tools/mash"):
    startt = time()
    check_call(
        f"{exe_path} triangle -p {threads} {sketch_file} > {dist_file}", shell=True
    )
    return time() - startt


def func_mash_search(
    query_sketch, ref_sketch, threads, dist_file, exe_path="./tools/mash"
):
    startt = time()
    check_call(
        f"{exe_path} dist {ref_sketch} {query_sketch} -p {threads} > {dist_file}",
        shell=True,
    )
    return time() - startt


def func_bindash_sketch(
    file_out, file_list, k, size, bits, threads, exe_path="./tools/bindash"
):
    start = time()
    check_call(
        f"{exe_path} sketch --minhashtype=2 --kmerlen={k} --sketchsize64={size//64} --bbits={bits} --nthreads={threads} --listfname={file_list} --outfname={file_out} 2> /tmp/tmp.dat",
        shell=True,
    )
    return time() - start


def func_bindash_dist(
    query_sketch, ref_sketch, out_file, threads, exe_path="./tools/bindash"
):
    start = time()
    check_call(
        f"{exe_path} dist {query_sketch} {ref_sketch} --nthreads={threads} --outfname={out_file} 2> /tmp/tmp.dat",
        shell=True,
    )
    return time() - start


def func_sourmash_sketch(file_out, file_list, k, scaled=1000, exe_path="sourmash"):
    starttm = time()
    check_call(
        f"{exe_path} sketch dna --from-file {file_list} -f --output {file_out} -p k={k},scaled={scaled},noabund",
        shell=True,
    )
    return time() - starttm


def func_sourmash_dist(
    query_sketch, ref_sketch, dist_fname, k, threads, scaled=1000, exe_path="sourmash"
):
    startt = time()
    if query_sketch is None:
        check_call(
            f"{exe_path} compare {ref_sketch} -p {threads} -k {k} --scaled {scaled} -f --ANI --csv {dist_fname}.csv --output {dist_fname}.mat",
            shell=True,
        )
    else:
        check_call(
            f"{exe_path} compare {query_sketch} {ref_sketch} -p {threads} -k {k} --scaled {scaled} -f --ANI --csv {dist_fname}.csv --output {dist_fname}.mat",
            shell=True,
        )
    return time() - startt


def func_dashing2_sketch(
    file_list, method, k, threads, size, exe_path="./tools/dashing2"
):
    start = time()
    check_call(
        f"{exe_path} sketch --{method} --cache -k {k} -S {size} -p {threads} -F {file_list}",
        shell=True,
    )
    # ./tools/dashing2 sketch --oph --cache -k 21 -S 1024 -F ./D1_10/fna_files.txt
    # ./tools/dashing2 sketch --probminhash --cache -k 21 -S 1024 -F ./D1_10/fna_files.txt
    # ./tools/dashing2 sketch --bagminhash --cache -k 21 -S 1024 -F ./D1_10/fna_files.txt
    return time() - start


def func_dashing2_dist(
    ref_file_list,
    query_file_list,
    out_file,
    method,
    k,
    size,
    threads,
    exe_path="./tools/dashing2",
):
    start = time()
    check_call(
        f"{exe_path} sketch --{method} --cache -k {k} -S {size} -p {threads} -F {ref_file_list} -Q {query_file_list} --cmpout {out_file}",
        shell=True,
    )
    return time() - start


def func_dashing2_search(
    query_file_list,
    ref_file_list,
    out_file,
    k,
    threads,
    size,
    exe_path="./tools/dashing2",
):
    start = time()
    check_call(
        f"{exe_path} cmp -k {k} -S {size} -p {threads} -F {ref_file_list} -Q {query_file_list} > {out_file}",
        shell=True,
    )
    return time() - start


def func_hd_sketch(
    fna_path,
    file_out,
    hash,
    k,
    S,
    D,
    threads,
    device="cpu",
    exe_path="hyper-gen",
):
    start = time()
    check_call(
        f"{exe_path} sketch -t {threads} -m {hash} -k {k} -s {S} -d {D} -D {device} -p {fna_path} -o {file_out} > tmp.dat",
        shell=True,
    )
    return time() - start


def func_hd_dist(
    query_sketch,
    ref_sketch,
    threads,
    dist_file,
    exe_path="hyper-gen",
):
    start = time()
    check_call(
        f"{exe_path} dist -t {threads} -r {ref_sketch} -q {query_sketch} -o {dist_file} > tmp.dat",
        shell=True,
    )
    return time() - start


def func_skani_sketch(
    file_out_path, file_list, c, m, threads, exe_path="./tools/skani"
):
    start = time()
    check_call(
        f"{exe_path} sketch  -t {threads} -c {c} -m {m} -l {file_list} -o {file_out_path}",
        shell=True,
    )
    return time() - start


def func_skani_dist(
    query_sketch, ref_sketch, dist_file, threads, exe_path="./tools/skani"
):
    start = time()
    check_call(
        f"{exe_path} dist -t {threads} -q {query_sketch} -r {ref_sketch} -o {dist_file} 2> /tmp/tmp.dat",
        shell=True,
    )
    return time() - start
