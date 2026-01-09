#!/usr/bin/env python3
"""
Contact-frequency postprocessing for GetContacts outputs.
Author: Yifang Chen

Step 1 (external tools; example)
-------------------------------
# Generate per-chain contact TSV
./getcontacts/get_dynamic_contacts.py \
  --topology ref.pdb --trajectory MD.xtc --itypes all \
  --output contacts_D_H.tsv --sele "chain D" --sele2 "chain H"

# Convert to per-residue-pair contact frequencies
./getcontacts/get_contact_frequencies.py \
  --input_files contacts_D_H.tsv --output_file resfrequencies_D_H.tsv

This script covers:
- Step 2: Convert 3-letter residue codes to 1-letter (optional).
- Step 3: Average 4 chain-pair TSVs within a repeat; filter + heatmap.
- Step 4: Average across repeats; filter + heatmap.

Example usage
-------------
# (Optional) convert residue names in a TSV
python contact_analysis.py convert resfrequencies_D_H.tsv converted_D_H.tsv

# Average 4 chain TSVs within a repeat + plot heatmap
python contact_analysis.py repeat \
  --inputs R1/converted_A_E_1.tsv R1/converted_B_F_1.tsv R1/converted_C_G_1.tsv R1/converted_D_H_1.tsv \
  --cutoff 0.2 --title RUN1 --outdir R1

# Average across repeats (using the per-repeat outputs from the previous command)
python contact_analysis.py average \
  --inputs R1/avg_contacts.tsv R2/avg_contacts.tsv \
  --cutoff 0.2 --title AVG --outdir .
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


RES3_TO_RES1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def convert_residue_name(res3: str) -> str:
    return RES3_TO_RES1.get(res3, res3)


def _convert_token(token: str) -> str:
    """
    Convert tokens like:
      'D:ASN:12' -> 'N12'   (drops chain prefix, converts ASN->N, keeps index)
      'H:GLY:5'  -> 'G5'
    If token doesn't match that pattern, return as-is.
    """
    token = re.sub(r"^[A-Za-z]:", "", token)  # drop leading "D:" / "H:" / etc.
    parts = token.split(":")
    if len(parts) == 2:
        # e.g., "ASN:12"
        res, idx = parts
        return f"{convert_residue_name(res)}{idx}"
    return token


def convert_tsv_resnames(input_file: Path, output_file: Path) -> None:
    with input_file.open("r", newline="") as fin, output_file.open("w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        for row in reader:
            writer.writerow([_convert_token(x) for x in row])


def read_contact_freq_tsv(path: Path) -> pd.DataFrame:
    """
    Expects GetContacts frequency TSV with 2 header lines.
    Columns: residue_1, residue_2, contact_frequency
    """
    return pd.read_csv(
        path, sep="\t", skiprows=2, header=None,
        names=["residue_1", "residue_2", "contact_frequency"],
    )


def average_contact_files(paths: Iterable[Path]) -> pd.DataFrame:
    paths = list(paths)
    if not paths:
        raise ValueError("No input files provided.")

    summed: defaultdict[Tuple[str, str], float] = defaultdict(float)
    for p in paths:
        df = read_contact_freq_tsv(p)
        for r1, r2, f in df[["residue_1", "residue_2", "contact_frequency"]].itertuples(index=False):
            summed[(r1, r2)] += float(f)

    n = len(paths)
    out = pd.DataFrame(
        [(r1, r2, freq / n) for (r1, r2), freq in summed.items()],
        columns=["residue_1", "residue_2", "avg_contact_frequency"],
    )

    # If duplicates exist, average them
    out = out.groupby(["residue_1", "residue_2"], as_index=False)["avg_contact_frequency"].mean()
    return out


def extract_residue_number(res_id: str) -> int:
    m = re.search(r"(\d+)", str(res_id))
    if not m:
        raise ValueError(f"Could not extract residue number from: {res_id}")
    return int(m.group(1))


def filter_and_sort(df: pd.DataFrame, cutoff: float) -> pd.DataFrame:
    df = df[df["avg_contact_frequency"] >= cutoff].copy()
    df["residue_1_num"] = df["residue_1"].map(extract_residue_number)
    df["residue_2_num"] = df["residue_2"].map(extract_residue_number)
    return df.sort_values(["residue_1_num", "residue_2_num"])


def pivot_for_heatmap(df: pd.DataFrame) -> pd.DataFrame:
    r1_order = df.drop_duplicates("residue_1").sort_values("residue_1_num")["residue_1"]
    r2_order = df.drop_duplicates("residue_2").sort_values("residue_2_num")["residue_2"]
    pivot = df.pivot(index="residue_1", columns="residue_2", values="avg_contact_frequency")
    return pivot.reindex(index=r1_order, columns=r2_order)


def plot_heatmap(
    pivot: pd.DataFrame,
    out_pdf: Path,
    title: str,
    xlabel: str = "Residue Index (LIP5)",
    ylabel: str = "Residue Index (AQP2)",
    vmin: float = 0.0,
    vmax: float = 1.0,
) -> None:
    plt.figure(figsize=(15, 10))
    ax = sns.heatmap(
        pivot,
        cmap="Blues",
        vmin=vmin,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        linecolor="gray",
        square=True,
        cbar_kws={"label": "Contact Frequency", "shrink": 0.5, "aspect": 10, "pad": 0.02},
    )
    cbar = ax.collections[0].colorbar
    cbar.ax.yaxis.label.set_size(14)

    plt.title(title, fontsize=16, weight="bold", pad=12)
    plt.xlabel(xlabel, fontsize=14, weight="bold")
    plt.ylabel(ylabel, fontsize=14, weight="bold")
    plt.xticks(fontsize=12, rotation=45, ha="right")
    plt.yticks(fontsize=12, rotation=45)
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.close()


def run_repeat(inputs: list[Path], outdir: Path, cutoff: float, title: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    avg_df = average_contact_files(inputs)
    avg_path = outdir / "avg_contacts.tsv"
    avg_df.to_csv(avg_path, sep="\t", index=False)

    filtered = filter_and_sort(avg_df.rename(columns={"avg_contact_frequency": "avg_contact_frequency"}), cutoff)
    filtered_path = outdir / "avg_contacts_cutoff_sorted.tsv"
    filtered.to_csv(filtered_path, sep="\t", index=False)

    pivot = pivot_for_heatmap(filtered)
    plot_heatmap(pivot, outdir / "contact_frequency_heatmap.pdf", title=title)


def run_average_repeat_outputs(inputs: list[Path], outdir: Path, cutoff: float, title: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    dfs = []
    for p in inputs:
        df = pd.read_csv(p, sep="\t")
        # Accept either "avg_contact_frequency" or "avg_contact_frequency" naming
        if "avg_contact_frequency" not in df.columns and "avg_contact_frequency" not in df.columns:
            # common expected column name from this script:
            if "avg_contact_frequency" not in df.columns and "avg_contact_frequency" not in df.columns:
                pass
        # normalize expected column name
        if "avg_contact_frequency" not in df.columns and "avg_contact_frequency" in df.columns:
            df = df.rename(columns={"avg_contact_frequency": "avg_contact_frequency"})
        if "avg_contact_frequency" in df.columns:
            df = df.rename(columns={"avg_contact_frequency": "avg_contact_frequency"})
        dfs.append(df[["residue_1", "residue_2", "avg_contact_frequency"]].copy())

    merged = dfs[0]
    for i, df in enumerate(dfs[1:], start=2):
        merged = merged.merge(df, on=["residue_1", "residue_2"], how="outer", suffixes=("", f"_r{i}"))

    freq_cols = [c for c in merged.columns if c.startswith("avg_contact_frequency")]
    merged["avg_contact_frequency"] = merged[freq_cols].mean(axis=1, skipna=True)

    out = merged[["residue_1", "residue_2", "avg_contact_frequency"]].copy()
    out = out.groupby(["residue_1", "residue_2"], as_index=False)["avg_contact_frequency"].mean()

    out_path = outdir / "avg_contacts_across_repeats.tsv"
    out.to_csv(out_path, sep="\t", index=False)

    filtered = filter_and_sort(out, cutoff)
    filtered_path = outdir / "avg_contacts_across_repeats_cutoff_sorted.tsv"
    filtered.to_csv(filtered_path, sep="\t", index=False)

    pivot = pivot_for_heatmap(filtered)
    plot_heatmap(pivot, outdir / "contact_frequency_avg_heatmap.pdf", title=title)


def main() -> None:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)

    p_conv = sub.add_parser("convert", help="Convert 3-letter residue codes to 1-letter in a TSV.")
    p_conv.add_argument("input", type=Path)
    p_conv.add_argument("output", type=Path)

    p_rep = sub.add_parser("repeat", help="Average multiple chain TSVs within a repeat and plot.")
    p_rep.add_argument("--inputs", nargs="+", type=Path, required=True)
    p_rep.add_argument("--outdir", type=Path, required=True)
    p_rep.add_argument("--cutoff", type=float, default=0.2)
    p_rep.add_argument("--title", type=str, default="REPEAT")

    p_avg = sub.add_parser("average", help="Average per-repeat outputs (avg_contacts.tsv) and plot.")
    p_avg.add_argument("--inputs", nargs="+", type=Path, required=True)
    p_avg.add_argument("--outdir", type=Path, required=True)
    p_avg.add_argument("--cutoff", type=float, default=0.2)
    p_avg.add_argument("--title", type=str, default="AVG")

    args = p.parse_args()

    if args.cmd == "convert":
        convert_tsv_resnames(args.input, args.output)
    elif args.cmd == "repeat":
        run_repeat(args.inputs, args.outdir, args.cutoff, args.title)
    elif args.cmd == "average":
        run_average_repeat_outputs(args.inputs, args.outdir, args.cutoff, args.title)


if __name__ == "__main__":
    main()
