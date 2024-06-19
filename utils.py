import os, glob, tqdm
import pandas as pd

import copy
import math
import numpy as np
from scipy import stats
import warnings


def get_all_files(path, extension):
    file_list = glob.glob(os.path.join(path, "*" + extension))
    return file_list


def parse_hd_params_from_filename(filename):
    params = filename.split("/")[-1].split(".")
    return {
        "alg": params[0],
        "dataset": params[1],
        "hash": params[3],
        "k": int(params[5]),
        "scaled": int(params[7]),
        "D": int(params[9]),
        "device": params[11],
    }


def parse_mash_params_from_filename(filename):
    params = filename.split("/")[-1].split(".")
    return {
        "alg": params[0],
        "dataset": params[2],
        "k": int(params[4]),
        "scaled": "-1",
        "D": int(params[6]),
    }


def parse_mash_search_file(mat_file):
    dist_dict = {}
    for line in open(mat_file, "r"):
        ref, query = line.split()[0:2]
        ref = ref.split("/")[-1]
        query = query.split("/")[-1]

        ani = (1 - float(line.split()[2])) * 100.0
        dist_dict[(ref, query)] = ani
    return dist_dict


def parse_bindash_params_from_filename(filename):
    params = filename.split("/")[-1].split(".")
    return {
        "alg": params[0],
        "dataset": params[2],
        "k": int(params[4]),
        "D": int(params[6]),
        "scaled": int(params[8]),
    }


def parse_bindash_dist_file(dist_file, ksize=21):
    dist_dict = {}
    for line in open(dist_file, "r"):
        ref, query = line.split()[0:2]
        ref = ref.split("/")[-1]
        query = query.split("/")[-1]

        J = line.split()[-1]
        J = float(J.split("/")[0]) / float(J.split("/")[1])
        if J > 1e-3:
            ANI = (
                max(
                    min(1 + 1 / float(ksize) * math.log(2 * J / (1 + J)), 1.0),
                    0.0,
                )
                * 100
            )

            dist_dict[(ref, query)] = ANI
    return dist_dict


def parse_matrix(mat_file, cutoff=0):
    fastani_add = 0.0
    count = 0
    dist_dict = dict()
    ref_vec = []
    dist_mat = []
    ani_cumulative = 0
    ani_counts = 0
    labels = []
    for line in open(mat_file, "r"):
        if count == 0:
            if "mash" in mat_file or "MASH" in mat_file:
                length = int(line.rstrip())
            elif len(line.split("\t")) > 2:
                length = len(line.split("\t"))
                labels = line.split("\t")
            else:
                length = int(line.split("\t")[-1].rstrip())
            dist_mat = np.zeros((length, length))
        elif count == 0 and "anim" in mat_file:
            length = len(line.split())
            dist_mat = np.zeros((length, length))
        if count != 0:
            spl = line.split()
            if "sour" in mat_file:
                ref = labels[count - 1].rstrip().split("/")[-1]
                endpoints = range(0, count - 1)
            else:
                ref = spl[0]
                ref = ref.split("/")[-1]
                endpoints = range(1, count)
            # if ".fa" not in ref:
            #     ref += ".fa"
            ref_vec.append(ref)
            for i in endpoints:
                try:
                    if "fastani" in mat_file:
                        ani = min(100, float(spl[i]) + fastani_add)
                    elif "mash" in mat_file or "MASH" in mat_file:
                        ani = (1 - float(spl[i])) * 100
                    else:
                        if float(spl[i]) <= 1:
                            ani = float(spl[i]) * 100
                        else:
                            ani = float(spl[i])
                    if "sour" in mat_file:
                        dist_mat[count - 1][i] = ani
                    else:
                        dist_mat[count - 1][i - 1] = ani
                    ani_cumulative += ani
                    ani_counts += 1
                except:
                    if "sour" in mat_file:
                        dist_mat[count - 1][i] = 0
                    else:
                        dist_mat[count - 1][i - 1] = 0

        count += 1
    for i in range(len(ref_vec)):
        for j in range(i):
            dist_dict[(ref_vec[i], ref_vec[j])] = dist_mat[i][j]
    if ani_cumulative / ani_counts > cutoff:
        return dist_dict
    else:
        return dict()


def parse_anim_excel(excel_file):
    df = pd.read_excel(excel_file).to_numpy()
    keys = [i[0] for i in df]

    dict_ani = {}
    for i, q in enumerate(keys):
        for j, r in enumerate(keys[i + 1 :]):
            dict_ani[(q + ".fna", r + ".fna")] = float(df[i, i + j + 2]) * 100
    return dict_ani


def parse_anim_csv(file):
    try:
        df = pd.read_csv(file, index_col=False)
        ani_dict = {
            (
                os.path.basename(df.loc[i, "query"]) + ".fna",
                df.loc[i, "reference"] + ".fna",
            ): float(df.loc[i, "ani"])
            * 100
            for i in range(len(df))
        }
    except:
        ani_dict = None
    return ani_dict


def convert_ani_to_jaccard(inp_dict, ksize):
    new_dict = {}
    for i in inp_dict:
        ani = inp_dict[i]
        new_dict[i] = 1 / (2 / math.exp((ani / 100 - 1) * ksize) - 1)
    return new_dict


def load_dash2_keys(key_file):
    # Read file
    with open(key_file, "r") as f:
        d = f.readlines()

    keys = []
    for i in range(len(d) // 2):
        q = d[2 * i][:-1] + ".fna"
        r = d[2 * i + 1][:-1] + ".fna"
        keys.append((q, r))
    return keys


def parse_dist_file(dist_file):
    dist_dict = dict()
    for line in open(dist_file, "r"):
        if "ANI" in line:
            continue
        spl = line.split()
        ref = spl[0]
        ref = ref.split("/")[-1]
        query = spl[1]
        query = query.split("/")[-1]
        if "fast" in dist_file or "anim" in dist_file:
            dist_dict[(query, ref)] = float(spl[2])
        elif "mash" in dist_file:
            dist_dict[(query, ref)] = (1 - float(spl[2])) * 100
        elif "aniu" in dist_file:
            dist_dict[(query, ref)] = float(spl[3])
        else:
            dist_dict[(query, ref)] = float(spl[2])

    return dist_dict


def parse_mash_dist_file(dist_file):
    dist_dict = parse_matrix(dist_file)
    return dist_dict


def parse_sourmash_ani_file(file):
    df = pd.read_csv(file)

    ani_dict = {}
    for i, fna in enumerate(df.columns):
        ref = os.path.basename(fna)

        for j in range(i):
            query = os.path.basename(df.columns[j])
            ani_dict[(query, ref)] = float(df.iloc[i, j]) * 100
    return ani_dict


def parse_dash2_params_from_filename(filename):
    params = filename.split("/")[-1].split(".")
    return {
        "alg": params[0],
        "hash": params[1],
        "dataset": params[2],
        "k": int(params[4]),
        "scaled": "-1",
        "D": int(params[6]),
    }


def parse_d2_dist_file(dist_file, ksize=0, if_ani=True):
    dist_dict = dict()
    f_read = open(dist_file, "r")
    dist_type = None
    for i, line in enumerate(f_read):
        if "#Dashing2 Symmetric pairwise Output" in line:
            dist_type = "pairwise"
            continue
        if "#Dashing2 Panel (Query/Refernce) Output" in line:
            dist_type = "search"
            continue

        if "#Sources" in line:
            keys = line.split()[1:]
            if dist_type == "search":
                ref_keys = [i.split("/")[-1] for i in keys[:-1]]
                query_key = keys[-1].split("/")[-1]

        if dist_type == "pairwise":
            if "#" in line or i == len(keys) + 2:
                continue

            spl = line.split()
            main_ref = spl[0].split("/")[-1]

            for j, J in enumerate(spl[i - 1 :]):
                other_ref = keys[i + j - 2].split("/")[-1]

                J = float(J)
                if J > 1e-3:
                    dist_dict[(main_ref, other_ref)] = (
                        (
                            max(
                                min(
                                    1 + 1 / float(ksize) * math.log(2 * J / (1 + J)),
                                    1.0,
                                ),
                                0.0,
                            )
                            * 100
                        )
                        if if_ani
                        else J
                    )

                # print(main_ref, other_ref, J, dist_dict[(main_ref, other_ref)])
        else:
            if "#" in line:
                continue
            spl = line.split()
            ref = spl[0].split("/")[-1]
            query = query_key

            J = float(spl[1])
            if J > 1e-3:
                dist_dict[(ref, query)] = (
                    (
                        max(
                            min(1 + 1 / float(ksize) * math.log(2 * J / (1 + J)), 1.0),
                            0.0,
                        )
                        * 100
                    )
                    if if_ani
                    else J
                )

    return dist_dict


def parse_skani_dist_file(dist_file):
    dist_dict = dict()
    for line in open(dist_file, "r"):
        spl = line.split()
        ref = spl[0].split("/")[-1]
        query = spl[1].split("/")[-1]
        if "ANI" in line:
            continue
        dist_dict[(ref, query)] = float(spl[2])

    return dist_dict


def parse_skani_params_from_filename(filename):
    params = filename.split("/")[-1].split(".")
    return {
        "alg": params[0],
        "dataset": params[2],
        "k": "-1",
        "scaled": int(params[6]),
        "D": int(params[4]),
    }


def report_error(ani_gt_est):
    ani_gt_est = np.array(ani_gt_est)
    mdive = np.mean((ani_gt_est[:, 1] - ani_gt_est[:, 0]) / np.max(ani_gt_est, axis=1))
    mae = np.mean(np.abs(ani_gt_est[:, 1] - ani_gt_est[:, 0]))
    rmse = np.sqrt(np.mean(np.square(ani_gt_est[:, 1] - ani_gt_est[:, 0])))
    mpe = np.mean((ani_gt_est[:, 1] - ani_gt_est[:, 0]) / ani_gt_est[:, 0]) * 100
    mpae = np.mean(np.abs(ani_gt_est[:, 1] - ani_gt_est[:, 0]) / ani_gt_est[:, 0]) * 100
    try:
        pearson_corr = stats.pearsonr(ani_gt_est[:, 0], ani_gt_est[:, 1])
        spearman_corr = stats.spearmanr(ani_gt_est[:, 0], ani_gt_est[:, 1])
        # print(
        #     "DivE = {:.4f}\tMAE = {:.4f}\tRMSE = {:.4f}\tMPE = {:.4f}\tMPAE = {:.4f}\nPearson = {} \nSpearman = {} for {} data points\n".format(
        #         mdive,
        #         mae,
        #         rmse,
        #         mpe,
        #         mpae,
        #         pearson_corr,
        #         spearman_corr,
        #         len(ani_gt_est),
        #     )
        # )

        return {
            "MAE": mae,
            "RMSE": rmse,
            "MPAE": mpae,
            "pearson": pearson_corr.statistic,
        }
    except:
        return {}


def compare_ani(ani_threshold, dict_ani_gt, dict_ani_cmp, desc=None, verbose=False):
    if desc:
        print(desc)

    processed_keys = []
    for k in list(dict_ani_cmp):
        k_reverse = (k[1], k[0])
        if (k_reverse in dict_ani_cmp) and (k_reverse not in processed_keys):
            processed_keys.append(k)
            dict_ani_cmp.pop(k_reverse)

    keys = dict_ani_cmp.keys()

    dict_ani_verbose = {"ref": [], "query": [], "ani_gt": [], "ani_est": []}
    ani_gt_est = []
    for k in keys:
        k_reverse = (k[1], k[0])

        if k in dict_ani_gt:
            ani_gt = dict_ani_gt[k]
        elif k_reverse in dict_ani_gt:
            ani_gt = dict_ani_gt[k_reverse]
        else:
            continue

        if (k in dict_ani_cmp) or (k_reverse in dict_ani_cmp):
            is_reverse = k_reverse in dict_ani_cmp

            ani_ = dict_ani_cmp[k_reverse] if is_reverse else dict_ani_cmp[k]

            if verbose:
                dict_ani_verbose["query"].append(k_reverse[0] if is_reverse else k[0])
                dict_ani_verbose["ref"].append(k_reverse[1] if is_reverse else k[1])
                dict_ani_verbose["ani_gt"].append(ani_gt)
                dict_ani_verbose["ani_est"].append(ani_)

            if ani_ >= ani_threshold and ani_gt >= ani_threshold:
                if ani_gt - ani_ > 10.0:
                    pass
                    # warnings.warn(
                    #     "{}: GT-Est ANI difference {} {} is greater than 10.0 for {} and {}".format(
                    #         desc, ani_gt, ani_, k[0], k[1]
                    #     )
                    # )
                else:
                    ani_gt_est.append([ani_gt, ani_])

    # print(desc, ani_gt_est)
    err = report_error(ani_gt_est)
    if verbose:
        return err, dict_ani_verbose
    else:
        return err, None


def compare_Jaccard(dict_J_gt, dict_J_cmp):
    keys = dict_J_gt.keys()

    jaccard_gt_est = []
    # dict_jaccard_verbose = {"ref": [], "query": [], "jaccard_gt": [], "jaccard_est": []}
    for k in keys:
        k_reverse = (k[1], k[0])

        if k in dict_J_gt:
            J_gt = dict_J_gt[k]
        elif k_reverse in dict_J_gt:
            J_gt = dict_J_gt[k_reverse]
        else:
            continue

        if (k in dict_J_cmp) or (k_reverse in dict_J_cmp):
            J_est = dict_J_cmp[k] if k in dict_J_cmp else dict_J_cmp[k_reverse]
            jaccard_gt_est.append([J_gt, J_est])

            # dict_jaccard_verbose["query"].append(k_reverse[0] if k_reverse else k[0])
            # dict_jaccard_verbose["ref"].append(k_reverse[1] if k_reverse else k[1])
            # dict_jaccard_verbose["ani_gt"].append(ani_gt)
            # dict_jaccard_verbose["ani_est"].append(ani_)

    err = report_error(jaccard_gt_est)
    return err
    # return dict_jaccard_verbose


def ani_to_jaccard(dict_ani, ksize):
    dict_jaccard = copy.copy(dict_ani)

    for i in dict_jaccard:
        ani = dict_jaccard[i]
        J = 1 / (2 / math.exp((ani / 100 - 1) * ksize) - 1)
        dict_jaccard[i] = J

    return dict_jaccard


def load_exp_file(file, parms):
    with open(file, "r") as f:
        lines = f.readlines()
        for i in lines[3:]:
            spl = i.split()
            t_sketch = float(spl[-2])
            t_search = float(spl[-1])
            if parms["alg"] == "HD":
                if (
                    parms["D"] == int(spl[3])
                    and parms["scaled"] == int(spl[4])
                    and parms["device"] == spl[6]
                ):
                    return {
                        "t_sketch": t_sketch,
                        "t_search": t_search,
                        "device": spl[-3],
                    }
            elif parms["alg"] in ["Dashing2", "bindash"]:
                if parms["D"] == int(spl[3]):
                    return {"t_sketch": t_sketch, "t_search": t_search}
            elif spl[0] in ["skani", "fastANI"]:
                return {"t_sketch": t_sketch, "t_search": t_search}
            elif parms["D"] == int(spl[2]):
                return {"t_sketch": t_sketch, "t_search": t_search}
