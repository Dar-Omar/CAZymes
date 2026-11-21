#!/usr/bin/env python3
"""Vertical bar plot of the number of CAZymes per substrate group.

A CAZyme is counted once in every group it belongs to.
A CAZyme with two substrates in the same row that belong to same group is counted once there.

Usage:
    python3 plots/plot_groups_total.py

Defaults:
    --input  dbcan_1071/overview_2_substrates_group.tsv
    --out    figures/groups_total.png
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter, MaxNLocator

PROJECT = Path(__file__).resolve().parents[2]


def project_path(p):
    """Relative paths are taken from the project directory."""
    p = Path(p).expanduser()
    return p if p.is_absolute() else PROJECT / p


def fmt(v, _):
    v = int(v)
    if v >= 1_000_000: return f"{v/1e6:g}M"
    if v >= 1_000: return f"{v/1e3:g}k"
    return str(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default="dbcan_1071/overview_2_substrates_group.tsv")
    ap.add_argument("--out", default="figures/groups_total.png")
    args = ap.parse_args()

    inp, out = project_path(args.input), project_path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # read only the group column; one row per CAZyme
    groups = pd.read_csv(inp, sep="\t", usecols=["Substrate_group"], dtype=str)["Substrate_group"]
    groups = groups[groups.notna() & (groups != "-")]
    tot = groups.str.split(";").explode().str.strip().value_counts()   # largest first

    if "UNMAPPED" in tot.index:
        print(f"WARNING: {tot['UNMAPPED']:,} CAZymes have an UNMAPPED group", file=sys.stderr)
    # put "Other" at the far right regardless of size
    if "Other" in tot.index:
        tot = pd.concat([tot.drop("Other"), tot[["Other"]]])

    fig, ax = plt.subplots(figsize=(0.8 * len(tot) + 1.5, 6))
    x = range(len(tot))
    ax.bar(x, tot.values, color="#9ecae1", width=0.7)
    ax.set_xticks(list(x))
    ax.set_xticklabels(tot.index, fontsize=10, rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_xlim(-0.6, len(tot) - 0.4)
    ymax = tot.max()
    ax.set_ylim(0, ymax * 1.08)
    for xx, v in zip(x, tot.values):
        ax.text(xx, v + ymax * 0.01, f"{v:,}", ha="center", va="bottom", fontsize=7.5)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    ax.yaxis.set_major_formatter(FuncFormatter(fmt))
    ax.set_ylabel("number of CAZymes", fontsize=11)
    ax.set_xlabel("CAZyme group", fontsize=11)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")


if __name__ == "__main__":
    main()
