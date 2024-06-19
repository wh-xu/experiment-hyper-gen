from argparse import ArgumentParser as AP
from sketchingfns import *
import sys
import os

acceptable_types = (
    "anim",
    "fastani",
    "bindash",
    "mash",
    "sourmash",
    "skani",
    "dashing2",
    "hd",
)


def main():
    ap = AP()
    ap.add_argument("--exp_type", choices=["search", "pw"], default="pw")
    ap.add_argument("--nrepeat", type=int, default=1)
    ap.add_argument("--nthreads", "-p", action="append", type=int)
    ap.add_argument("-k", nargs="+", type=int)
    ap.add_argument("--method", nargs="+", type=str)
    ap.add_argument("--hd_scaled", nargs="+", type=int)
    ap.add_argument("--hd_dim", nargs="+", type=int)
    ap.add_argument("--device", nargs="+", type=str)
    ap.add_argument("--sketchsize", "-s", nargs="+", type=int)
    ap.add_argument(
        "--only_type",
        help="If set, only performs benchmark for a particular type.",
        type=str,
    )
    ap.add_argument("-query_list", required=True, type=str)
    ap.add_argument("-ref_list", required=True, type=str)
    ap.add_argument("-suffix", required=True, type=str, help="Descriptive file suffix")

    args = ap.parse_args()
    exp_type = args.exp_type
    query_file_list, ref_file_list = args.query_list, args.ref_list
    assert os.path.isfile(
        query_file_list
    ), f"fnames argument ({query_file_list}) does not exist"
    assert os.path.isfile(
        ref_file_list
    ), f"fnames argument ({ref_file_list}) does not exist"

    if not args.sketchsize:
        raise ValueError("Required: at least one value for sketchsize")
    if not args.k:
        raise ValueError("Required: at least one value for k")

    sszes = args.sketchsize
    if any(x <= 0 for x in sszes):
        raise ValueError("sketchsize must be > 1 and a power of 2")

    method = args.method
    hd_dims = args.hd_dim
    scaled = args.hd_scaled
    device = args.device

    print(f"##Options: {args}")
    print(f"##Results for {query_file_list} vs. {ref_file_list}")
    print(
        "#Method\tk\tNumReg\tRegisterSize\tNumThreads\tSketchTime\tDistanceTime",
        flush=True,
    )

    # read a txt file with multiple lines and store into a list
    with open(query_file_list, "r") as f:
        query_files = f.read().splitlines()

    with open(ref_file_list, "r") as f:
        ref_files = f.read().splitlines()

    rstr = args.suffix
    for nt in args.nthreads:
        if nt < 0:
            from multiprocessing import cpu_count

            nt = cpu_count()

        # 1. ANIm
        if "anim" == args.only_type:
            dir_out = f"results_anim/anim.{rstr}"
            trun = repeat_x(
                func_run_anim,
                args.nrepeat,
                query_files=query_files,
                ref_files=ref_files,
                path_out=dir_out,
                threads=nt,
            )
            print(f"ANIm\t-\t-\t-\t{nt}\t-\t{trun}", flush=True)

        if "skani" == args.only_type:
            path_out = "results_skani"

            c, m = 70, 1000
            if exp_type == "search":
                _ = repeat_x(
                    func_skani_sketch,
                    args.nrepeat,
                    file_list=query_file_list,
                    file_out_path=f"{path_out}/{rstr}/query",
                    c=c,
                    m=m,
                    threads=nt,
                )

                tsketch = repeat_x(
                    func_skani_sketch,
                    args.nrepeat,
                    file_list=ref_file_list,
                    file_out_path=f"{path_out}/{rstr}/ref",
                    c=c,
                    m=m,
                    threads=nt,
                )

                query_sketch = f"./{path_out}/{rstr}/query/*.sketch"
                ref_sketch = f"./{path_out}/{rstr}/ref/*.sketch"
                dist_file = f"./{path_out}/skani.dist.{rstr}.c.{c}.m.{m}.dist"
                tdist = repeat_x(
                    func_skani_dist,
                    args.nrepeat,
                    query_sketch=query_sketch,
                    ref_sketch=ref_sketch,
                    dist_file=dist_file,
                    threads=nt,
                )
            else:
                tsketch = repeat_x(
                    func_skani_sketch,
                    args.nrepeat,
                    file_out_path=f"{path_out}/{rstr}",
                    file_list=ref_file_list,
                    c=c,
                    m=m,
                    threads=nt,
                )

                sketch_file = f"./{path_out}/{rstr}/*.sketch"
                dist_file = f"./{path_out}/skani.dist.{rstr}.c.{c}.m.{m}.dist"
                tdist = repeat_x(
                    func_skani_dist,
                    args.nrepeat,
                    query_sketch=sketch_file,
                    ref_sketch=sketch_file,
                    dist_file=dist_file,
                    threads=nt,
                )

            print(f"skani\t-\t{c}\t{m}\t{nt}\t{tsketch}\t{tdist}", flush=True)

        # 2. fastANI
        if "fastani" == args.only_type:
            path_out = "results_fastani"
            os.makedirs(path_out, exist_ok=True)

            # Split large refs into batches
            batch_size = 2000
            with open(ref_file_list, "r") as f:
                file_read = f.readlines()
                num_fna_ref = len(file_read)

            if num_fna_ref > batch_size:
                num_batch = (num_fna_ref + batch_size) // batch_size

                for i in range(num_batch):
                    _ref_file_list = ref_file_list[:-4] + f"_{i}.txt"
                    with open(_ref_file_list, "w") as f:
                        print(f"Splitting files into {i}/{num_batch} sub-batches")
                        f.write(
                            "".join(file_read[i * batch_size : (i + 1) * batch_size])
                        )

                    dist_file = f"fastani.{rstr}_{i}.log"
                    trun = repeat_x(
                        func_run_fastani,
                        args.nrepeat,
                        query_files=query_file_list,
                        ref_files=_ref_file_list,
                        path_out=path_out,
                        file_out=dist_file,
                        threads=nt,
                    )
                    print(f"fastANI-batch-{i}\t-\t-\t-\t{nt}\t0\t{trun}", flush=True)
                check_call(
                    f"cat {path_out}/fastani.{rstr}_*.log > {path_out}/fastani.{rstr}.log",
                    shell=True,
                )
            else:
                dist_file = f"fastani.{rstr}.log"
                trun = repeat_x(
                    func_run_fastani,
                    args.nrepeat,
                    query_files=query_file_list,
                    ref_files=ref_file_list,
                    path_out=path_out,
                    file_out=dist_file,
                    threads=nt,
                )
                print(f"fastANI\t-\t-\t-\t{nt}\t0\t{trun}", flush=True)

        for k in args.k:
            # 4. Sourmash
            if "sourmash" == args.only_type:
                path_out = "results_sourmash"
                os.makedirs(path_out, exist_ok=True)

                scaled = 1000
                if exp_type == "search":
                    pass
                else:
                    sketch_file = (
                        f"./{path_out}/SOURMASH.sketch.{rstr}.k.{k}.scaled.{scaled}.sig"
                    )
                    tsketch = repeat_x(
                        func_sourmash_sketch,
                        args.nrepeat,
                        file_out=sketch_file,
                        file_list=ref_file_list,
                        k=k,
                        scaled=scaled,
                    )

                    os.makedirs(f"./{path_out}/{rstr}", exist_ok=True)
                    dist_file = f"./{path_out}/{rstr}/SOURMASH.dist.{rstr}.k.{k}.scaled.{scaled}"
                    tdist = repeat_x(
                        func_sourmash_dist,
                        args.nrepeat,
                        query_sketch=None,
                        ref_sketch=sketch_file,
                        dist_fname=dist_file,
                        k=k,
                        threads=nt,
                        scaled=scaled,
                    )
                    print(
                        f"Sourmash\t{k}\t{scaled}\t8\t{nt}\t{tsketch}\t{tdist}",
                        flush=True,
                    )

                # sourmash search query.sig [ list of signatures or SBTs ]

            # 4. HD
            if args.only_type == "hd":
                path_out = "results_hd_rs"
                os.makedirs(path_out, exist_ok=True)
                os.makedirs(f"{path_out}/{rstr}", exist_ok=True)

                for dev in device:
                    for S in scaled:
                        for D in hd_dims:
                            for hash in method:
                                if exp_type == "search":
                                    query_path = (
                                        f"../hyper-gen-rust/fna_query/{rstr[:-4]}"
                                    )
                                    query_sketch_file = f"{path_out}/{rstr}/HD.query.h.{hash}.k.{k}.s.{S}.D.{D}.device.{dev}.sketch"
                                    _ = repeat_x(
                                        func_hd_sketch,
                                        args.nrepeat,
                                        fna_path=query_path,
                                        file_out=query_sketch_file,
                                        hash=hash,
                                        k=k,
                                        S=S,
                                        D=D,
                                        device=dev,
                                        threads=nt,
                                    )

                                    ref_sketch_file = f"./{path_out}/{rstr}/HD.{rstr}.h.{hash}.k.{k}.s.{S}.D.{D}.device.{dev}.sketch"
                                    ref_path = os.path.dirname(ref_file_list)
                                    tsketch = repeat_x(
                                        func_hd_sketch,
                                        args.nrepeat,
                                        fna_path=ref_path,
                                        file_out=ref_sketch_file,
                                        hash=hash,
                                        k=k,
                                        S=S,
                                        D=D,
                                        device=dev,
                                        threads=nt,
                                    )

                                    dstfile = f"./{path_out}/{rstr}/HD.{rstr}.h.{hash}.k.{k}.s.{S}.D.{D}.device.{dev}.dist"
                                    tdist = repeat_x(
                                        func_hd_dist,
                                        args.nrepeat,
                                        query_sketch=query_sketch_file,
                                        ref_sketch=ref_sketch_file,
                                        threads=nt,
                                        dist_file=dstfile,
                                    )
                                else:
                                    ref_sketch_file = f"./{path_out}/{rstr}/HD.{rstr}.h.{hash}.k.{k}.s.{S}.D.{D}.device.{dev}.sketch"
                                    ref_path = os.path.dirname(ref_file_list)
                                    tsketch = repeat_x(
                                        func_hd_sketch,
                                        args.nrepeat,
                                        fna_path=ref_path,
                                        file_out=ref_sketch_file,
                                        hash=hash,
                                        k=k,
                                        S=S,
                                        D=D,
                                        device=dev,
                                        threads=nt,
                                    )

                                    dstfile = f"./{path_out}/{rstr}/HD.{rstr}.h.{hash}.k.{k}.s.{S}.D.{D}.device.{dev}.dist"
                                    tdist = repeat_x(
                                        func_hd_dist,
                                        args.nrepeat,
                                        query_sketch=ref_sketch_file,
                                        ref_sketch=ref_sketch_file,
                                        threads=nt,
                                        dist_file=dstfile,
                                    )

                                print(
                                    f"HD\t{hash}\t{k}\t{D}\t{S}\t{nt}\t{dev}\t{tsketch:.4f}\t{tdist:.4f}",
                                    flush=True,
                                )

            for ssz in sszes:
                # 3. MASH
                if "mash" == args.only_type:
                    path_out = "results_mash"
                    os.makedirs(path_out, exist_ok=True)

                    if exp_type == "search":
                        query_sketch_file = (
                            f"./{path_out}/MASH.querysketch.{rstr}.k.{k}.S.{ssz}.msh"
                        )
                        _ = repeat_x(
                            func_mash_sketch,
                            args.nrepeat,
                            file_out=query_sketch_file,
                            file_list=query_file_list,
                            k=k,
                            size=ssz,
                            threads=nt,
                        )

                        ref_sketch_file = (
                            f"./{path_out}/MASH.sketch.{rstr}.k.{k}.S.{ssz}.msh"
                        )
                        # if not os.path.isfile(ref_sketch_file):
                        tsketch = repeat_x(
                            func_mash_sketch,
                            args.nrepeat,
                            file_out=ref_sketch_file,
                            file_list=ref_file_list,
                            k=k,
                            size=ssz,
                            threads=nt,
                        )
                        # else:
                        # tsketch = -1

                        mdstfile = f"./{path_out}/MASH.dist.{rstr}.k.{k}.S.{ssz}.phylip"
                        tdist = repeat_x(
                            func_mash_search,
                            args.nrepeat,
                            query_sketch=query_sketch_file,
                            ref_sketch=ref_sketch_file,
                            threads=nt,
                            dist_file=mdstfile,
                        )
                    else:
                        sketch_file = (
                            f"./{path_out}/MASH.sketch.{rstr}.k.{k}.S.{ssz}.msh"
                        )
                        tsketch = repeat_x(
                            func_mash_sketch,
                            args.nrepeat,
                            file_out=sketch_file,
                            file_list=ref_file_list,
                            k=k,
                            size=ssz,
                            threads=nt,
                        )

                        mdstfile = f"./{path_out}/MASH.dist.{rstr}.k.{k}.S.{ssz}.phylip"
                        tdist = repeat_x(
                            func_mash_dist,
                            args.nrepeat,
                            sketch_file=sketch_file,
                            threads=nt,
                            dist_file=mdstfile,
                        )

                    print(f"Mash\t{k}\t{ssz}\t8\t{nt}\t{tsketch}\t{tdist}", flush=True)

                # 5. bindash
                if "bindash" == args.only_type:
                    path_out = "results_bindash"
                    os.makedirs(path_out, exist_ok=True)

                    bits = 14
                    dist_file = f"./{path_out}/bindash.dist.{rstr}.k.{k}.S.{ssz}.bits.{bits}.dist"
                    if exp_type == "search":
                        ref_sketch_file = f"./{path_out}/bindash.sketch.{rstr}.k.{k}.s.{ssz}.bits.{bits}.sketch"
                        tsketch = repeat_x(
                            func_bindash_sketch,
                            args.nrepeat,
                            file_out=ref_sketch_file,
                            file_list=ref_file_list,
                            k=k,
                            bits=bits,
                            size=ssz,
                            threads=nt,
                        )

                        query_sketch_file = f"./{path_out}/bindash.querysketch.{rstr}.k.{k}.s.{ssz}.bits.{bits}.sketch"
                        _ = repeat_x(
                            func_bindash_sketch,
                            args.nrepeat,
                            file_out=query_sketch_file,
                            file_list=query_file_list,
                            k=k,
                            bits=bits,
                            size=ssz,
                            threads=nt,
                        )

                        tdist = repeat_x(
                            func_bindash_dist,
                            args.nrepeat,
                            query_sketch=query_sketch_file,
                            ref_sketch=ref_sketch_file,
                            out_file=dist_file,
                            threads=nt,
                        )
                    else:
                        ref_sketch_file = f"./{path_out}/bindash.sketch.{rstr}.k.{k}.s.{ssz}.bits.{bits}.sketch"
                        tsketch = repeat_x(
                            func_bindash_sketch,
                            args.nrepeat,
                            file_out=ref_sketch_file,
                            file_list=ref_file_list,
                            k=k,
                            bits=bits,
                            size=ssz,
                            threads=nt,
                        )

                        tdist = repeat_x(
                            func_bindash_dist,
                            args.nrepeat,
                            ref_sketch=ref_sketch_file,
                            query_sketch=ref_sketch_file,
                            out_file=dist_file,
                            threads=nt,
                        )
                    print(
                        f"bindash\t{k}\t{bits}\t{ssz}\t\t{nt}\t{tsketch}\t{tdist}",
                        flush=True,
                    )

                # 6. Dashing 2
                if "dashing2" == args.only_type:
                    path_out = "results_dashing2"
                    os.makedirs(path_out, exist_ok=True)

                    for hash in method:
                        dist_file = (
                            f"./{path_out}/Dashing2.{hash}.{rstr}.k.{k}.S.{ssz}.dist"
                        )
                        if exp_type == "search":
                            tsketch = repeat_x(
                                func_dashing2_sketch,
                                args.nrepeat,
                                file_list=ref_file_list,
                                method=hash,
                                k=k,
                                size=ssz,
                                threads=nt,
                            )

                            _ = repeat_x(
                                func_dashing2_sketch,
                                args.nrepeat,
                                file_list=query_file_list,
                                method=hash,
                                k=k,
                                size=ssz,
                                threads=nt,
                            )

                            tdist = repeat_x(
                                func_dashing2_dist,
                                args.nrepeat,
                                ref_file_list=ref_file_list,
                                query_file_list=query_file_list,
                                out_file=dist_file,
                                method=hash,
                                k=k,
                                size=ssz,
                                threads=nt,
                            )
                        else:
                            tsketch = repeat_x(
                                func_dashing2_sketch,
                                args.nrepeat,
                                file_list=ref_file_list,
                                method=hash,
                                k=k,
                                size=ssz,
                                threads=nt,
                            )

                            tdist = repeat_x(
                                func_dashing2_dist,
                                args.nrepeat,
                                ref_file_list=ref_file_list,
                                query_file_list=ref_file_list,
                                out_file=dist_file,
                                method=hash,
                                k=k,
                                size=ssz,
                                threads=nt,
                            )
                        print(
                            f"Dashing2\t{hash}\t{k}\t{ssz}\t\t{nt}\t{tsketch}\t{tdist}",
                            flush=True,
                        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
