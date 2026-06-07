"""
Generators for figure6 (spatial), figure8 (sexual conservation), figure9
(enrichment) — previously committed as PNGs with no source script. All three
are now reproduced from network_enriched.csv + results/derived_stats.json
(+ BANC metadata for positions).

Run:  MCIS_DATA_DIR=/path/to/data python src/figures_extra.py
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mcis_paths import data_dir, repo_root
REPO = repo_root()
FIG = os.path.join(REPO, "figures")
DATA_DIR = data_dir()

CLASS_COLORS = {
    "descending": "#d62728", "ascending": "#1f77b4",
    "sensory_ascending": "#2ca02c", "sensory_descending": "#9467bd",
}


def _style(ax, title):
    ax.set_title(title, fontsize=11, fontweight="bold", pad=6)
    ax.spines[["top", "right"]].set_visible(False)


def figure6_spatial(enr, meta):
    pos = meta.set_index("root_626")["root_position_nm"].dropna()

    def xyz(ids):
        out = []
        for i in ids:
            if i in pos.index:
                try:
                    out.append([float(v) for v in str(pos[i]).split(",")])
                except ValueError:
                    pass
        return np.array(out) / 1000.0  # nm → µm

    circ = xyz(enr["BANC"].astype(str))
    bg = xyz(meta[meta["fafb_match"].notna() & meta["manc_match"].notna()]
             ["root_626"].astype(str))
    classes = enr["super_class"].tolist()
    cols = [CLASS_COLORS.get(c, "#888") for c in classes]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor("white")
    for ax, (a, b, name) in zip(axes, [(0, 1, "coronal (X-Y)"),
                                       (0, 2, "sagittal (X-Z)"),
                                       (1, 2, "axial (Y-Z)")]):
        if len(bg):
            ax.scatter(bg[:, a], bg[:, b], s=3, c="#ddd", alpha=0.4)
        if len(circ):
            ax.scatter(circ[:, a], circ[:, b], s=22, c=cols, edgecolors="k",
                       linewidths=0.3)
        ax.set_xlabel(f"{name} (µm)"); _style(ax, name)
        ax.set_aspect("equal", "datalim")
    handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=v, mec="k", label=k)
               for k, v in CLASS_COLORS.items()]
    axes[0].legend(handles=handles, fontsize=8, loc="best")
    fig.suptitle(f"Anatomical distribution of the {len(enr)}-neuron circuit "
                 "(grey = all matched neurons)", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "figure6_spatial.png"), dpi=150,
                bbox_inches="tight", facecolor="white")
    print("Saved figure6_spatial.png")


def figure8_sexual(enr):
    n = len(enr)
    iso = int((enr["sexually_dimorphic"] == "isomorphic").sum())
    dim = n - iso
    dimorphic = enr[enr["sexually_dimorphic"] != "isomorphic"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.patch.set_facecolor("white")
    axes[0].pie([iso, dim], labels=[f"isomorphic\n{iso}", f"dimorphic\n{dim}"],
                colors=["#2ca02c", "#d62728"], autopct="%1.1f%%",
                startangle=90, wedgeprops=dict(width=0.45))
    _style(axes[0], f"A  Sexual conservation ({100*iso/n:.1f}%)")

    by_cls = dimorphic["super_class"].value_counts()
    axes[1].bar(by_cls.index, by_cls.values,
                color=[CLASS_COLORS.get(c, "#888") for c in by_cls.index])
    axes[1].set_ylabel("dimorphic count")
    axes[1].tick_params(axis="x", rotation=30)
    _style(axes[1], "B  Dimorphism by class")

    nt = dimorphic["neurotransmitter_predicted"].value_counts()
    axes[2].bar(nt.index, nt.values, color="#9467bd")
    axes[2].set_ylabel("count"); axes[2].tick_params(axis="x", rotation=30)
    _style(axes[2], "C  NT of dimorphic neurons")
    fig.suptitle("Sexual conservation of the conserved circuit",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "figure8_sexual_conservation.png"), dpi=150,
                bbox_inches="tight", facecolor="white")
    print("Saved figure8_sexual_conservation.png")


def figure9_enrichment(enr, derived):
    enrich = derived["enrichment"]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.patch.set_facecolor("white")

    cls = list(enrich.keys())
    folds = [enrich[c]["fold"] for c in cls]
    axes[0].bar(cls, folds, color=[CLASS_COLORS.get(c, "#888") for c in cls])
    for i, c in enumerate(cls):
        axes[0].text(i, folds[i], f"{folds[i]:.1f}×", ha="center", va="bottom")
    axes[0].set_ylabel("fold enrichment vs FAFB")
    _style(axes[0], "A  Superclass enrichment")

    nt = enr["neurotransmitter_predicted"].value_counts()
    axes[1].bar(nt.index, nt.values, color="#1f77b4")
    axes[1].set_ylabel("count"); axes[1].tick_params(axis="x", rotation=30)
    _style(axes[1], "B  Neurotransmitter profile")

    hl = enr["hemilineage"].value_counts().head(6)
    axes[2].barh(hl.index[::-1], hl.values[::-1], color="#2ca02c")
    axes[2].set_xlabel("count")
    _style(axes[2], "C  Top hemilineages")
    fig.suptitle("Cell-type enrichment vs FAFB whole-brain background",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "figure9_enrichment.png"), dpi=150,
                bbox_inches="tight", facecolor="white")
    print("Saved figure9_enrichment.png")


def main():
    enr = pd.read_csv(os.path.join(REPO, "network_enriched.csv"),
                      dtype={"BANC": str})
    with open(os.path.join(REPO, "results", "derived_stats.json")) as f:
        derived = json.load(f)
    meta = pd.read_feather(DATA_DIR + "banc_meta.feather")
    meta["root_626"] = meta["root_626"].astype(str)
    meta = meta.drop_duplicates("root_626")
    figure6_spatial(enr, meta)
    figure8_sexual(enr)
    figure9_enrichment(enr, derived)


if __name__ == "__main__":
    main()
