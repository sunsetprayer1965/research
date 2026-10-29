#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import numpy as np
import pandas as pd


def pick_one(ds_dir: str, kind: str) -> str:
    patterns = [
        os.path.join(ds_dir, f"interactions_*_{kind}*.txt"),
        os.path.join(ds_dir, f"interactions_*_{kind}*.txt.txt"),
    ]
    matches = []
    for p in patterns:
        matches.extend(glob.glob(p))
    matches = sorted(set(matches))
    if not matches:
        return ""
    if kind == "intra":
        all_matches = [m for m in matches if "all" in os.path.basename(m).lower()]
        if all_matches:
            return all_matches[0]
    return matches[0]


def load_locus_assignment(bed_path: str) -> pd.DataFrame:
    locus_assignment = pd.read_csv(
        bed_path,
        delim_whitespace=True,
        names=["chromosome", "pair1", "pair2", "locus"],
    )

    locus_assignment["pair1"] = pd.to_numeric(locus_assignment["pair1"], errors="coerce").fillna(0).astype(int)
    locus_assignment["pair2"] = pd.to_numeric(locus_assignment["pair2"], errors="coerce").fillna(0).astype(int)

    locus_assignment["chromosome"] = (
        locus_assignment["chromosome"]
        .astype(str)
        .str.replace("chr", "", regex=False)
    )
    locus_assignment["chromosome"] = pd.to_numeric(locus_assignment["chromosome"], errors="coerce").astype(int)

    locus_assignment["locus"] = pd.to_numeric(locus_assignment["locus"], errors="coerce")

    # 防御：清理异常行
    locus_assignment = locus_assignment.dropna(subset=["chromosome", "pair1", "pair2", "locus"])
    locus_assignment = locus_assignment.astype({"chromosome": int, "pair1": int, "pair2": int})

    return locus_assignment


def load_interactions(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        names=["chr1", "locus1", "chr2", "locus2", "freq", "qvalue"],
        skiprows=1,
    )
    for c in ["chr1", "locus1", "chr2", "locus2", "freq"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["chr1", "locus1", "chr2", "locus2", "freq"])
    df["chr1"] = df["chr1"].astype(int)
    df["chr2"] = df["chr2"].astype(int)
    df["locus1"] = df["locus1"].astype(int)
    df["locus2"] = df["locus2"].astype(int)
    return df


def build_lookup(locus_assignment: pd.DataFrame):
    """
    Build fast lookup:
      lookup[chrom] -> numpy array mapping pair_index -> locus (float), NaN if missing
    """
    print("[locus] building lookup...", flush=True)

    lookup = {}

    for chrom, g in locus_assignment.groupby("chromosome", sort=True):
        max_pair = int(g["pair2"].max())
        arr = np.full(max_pair + 1, np.nan, dtype=float)

        # 填充区间：pair1..pair2 -> locus
        # 这里用 for 循环遍历 bed 行数（通常远小于 interactions 行数，所以很快）
        for p1, p2, loc in g[["pair1", "pair2", "locus"]].itertuples(index=False, name=None):
            if p2 >= 0:
                p1 = max(p1, 0)
                arr[p1 : p2 + 1] = loc

        lookup[int(chrom)] = arr

    print(f"[locus] lookup built for {len(lookup)} chromosomes", flush=True)
    return lookup


def map_with_lookup(df: pd.DataFrame, lookup, name: str) -> pd.DataFrame:
    """
    Vectorized mapping using lookup arrays.
    """
    print(f"[locus] mapping {name} with lookup...", flush=True)

    # chr1 mapping
    locus1_mapped = np.full(len(df), np.nan, dtype=float)
    locus2_mapped = np.full(len(df), np.nan, dtype=float)

    # 按 chr 分组，避免逐行查字典
    for chrom in sorted(set(df["chr1"].unique())):
        arr = lookup.get(int(chrom))
        if arr is None:
            continue
        idx = (df["chr1"].values == chrom)
        pos = df.loc[idx, "locus1"].values
        ok = (pos >= 0) & (pos < arr.shape[0])
        tmp = np.full(pos.shape[0], np.nan, dtype=float)
        tmp[ok] = arr[pos[ok]]
        locus1_mapped[idx] = tmp

    for chrom in sorted(set(df["chr2"].unique())):
        arr = lookup.get(int(chrom))
        if arr is None:
            continue
        idx = (df["chr2"].values == chrom)
        pos = df.loc[idx, "locus2"].values
        ok = (pos >= 0) & (pos < arr.shape[0])
        tmp = np.full(pos.shape[0], np.nan, dtype=float)
        tmp[ok] = arr[pos[ok]]
        locus2_mapped[idx] = tmp

    out = pd.DataFrame({
        "locus1": locus1_mapped,
        "locus2": locus2_mapped,
        "frequency": df["freq"].values.astype(float),
    })

    out = out.dropna(subset=["locus1", "locus2"])
    print(f"[locus] {name} mapped rows: {len(out)}", flush=True)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset_dir", required=True)
    parser.add_argument("--bed", default="duan.SC.10000.raw_sub.bed")
    parser.add_argument("--out", default="duan.SC.10000.raw_sub.matrix")
    args = parser.parse_args()

    print("[locus] START", flush=True)
    print(f"[locus] dataset_dir={args.dataset_dir}", flush=True)
    print(f"[locus] bed={args.bed}", flush=True)
    print(f"[locus] out={args.out}", flush=True)

    inter_file = pick_one(args.dataset_dir, "inter")
    intra_file = pick_one(args.dataset_dir, "intra")
    if not inter_file or not intra_file:
        raise FileNotFoundError(f"[locus] Missing inter/intra in {args.dataset_dir}")

    print(f"[locus] inter_file = {inter_file}", flush=True)
    print(f"[locus] intra_file = {intra_file}", flush=True)

    print("[locus] loading bed...", flush=True)
    locus_assignment = load_locus_assignment(args.bed)
    print(f"[locus] bed loaded rows={len(locus_assignment)}", flush=True)

    lookup = build_lookup(locus_assignment)

    print("[locus] loading inter...", flush=True)
    inter_df = load_interactions(inter_file)
    print(f"[locus] inter loaded rows={len(inter_df)}", flush=True)

    print("[locus] loading intra...", flush=True)
    intra_df = load_interactions(intra_file)
    print(f"[locus] intra loaded rows={len(intra_df)}", flush=True)

    mapped_inter = map_with_lookup(inter_df, lookup, "inter")
    mapped_intra = map_with_lookup(intra_df, lookup, "intra")

    combined = pd.concat([mapped_inter, mapped_intra], ignore_index=True)
    final = combined.groupby(["locus1", "locus2"], as_index=False)["frequency"].sum()

    final = final[final["locus1"] != final["locus2"]].copy()
    final["frequency"] = final["frequency"].astype(float)

    final.to_csv(args.out, sep="\t", index=False, header=False)

    print(f"[locus] wrote: {os.path.abspath(args.out)}", flush=True)
    print(f"[locus] rows: {len(final)}", flush=True)


if __name__ == "__main__":
    main()


# 1) 我是怎么把 locus.py 提速的（简单总结）

# 你原来的慢点在这里：

# 对 interactions 的每一行都 apply(axis=1)

# 每一行都去 locus_assignment DataFrame 里做一次过滤：
# (chrom==...) & (pair1<=idx) & (pair2>=idx)
# 这等于“Python 循环 × DataFrame筛选”，遇到你这种百万级行数会非常慢。

# 我做的提速核心是：

# 用 bed 的区间信息先构建一个查表结构（lookup）：
# lookup[chrom][pair_index] = locus（一个 numpy 数组）

# mapping 时不再逐行筛 DataFrame，而是向量化查表：
# locus1_mapped = lookup[chr1][locus1]
# locus2_mapped = lookup[chr2][locus2]

# 这样 mapping 变成 O(N) 的数组索引，通常能快一个数量级以上（你刚刚 log 也验证了）。