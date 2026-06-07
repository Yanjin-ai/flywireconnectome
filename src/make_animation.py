"""
Narrative 3D animation of the conserved circuit (Line C completion).
====================================================================

Renders a rotating 3D view of the conserved circuit in BANC anatomical space,
neurons coloured by their conservation z-score and the conserved edges drawn in
gold. Saves an animated GIF (matplotlib + Pillow; no Blender/manim needed).

Run:  MCIS_DATA_DIR=/path/to/data python src/make_animation.py
Output: figures/circuit_3d_conservation.gif
"""
import os
import sys

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, figures_dir, repo_root, results_dir  # noqa


def main():
    repo = repo_root()
    dd = data_dir()
    enr = pd.read_csv(os.path.join(repo, "network_enriched.csv"),
                      dtype={"BANC": str})
    cons_p = os.path.join(results_dir(), "neuron_conservation.csv")
    if os.path.exists(cons_p):
        cons = pd.read_csv(cons_p, dtype={"BANC": str})
        enr = enr.merge(cons[["BANC", "conservation_z"]], on="BANC", how="left")
    else:
        enr["conservation_z"] = 0.0

    meta = pd.read_feather(dd + "banc_meta.feather")
    meta["root_626"] = meta["root_626"].astype(str)
    pos = (meta.drop_duplicates("root_626")
           .set_index("root_626")["root_position_nm"].dropna())

    P, z, banc_ok = [], [], []
    for _, r in enr.iterrows():
        if r["BANC"] in pos.index:
            try:
                P.append([float(x) for x in str(pos[r["BANC"]]).split(",")])
                z.append(r["conservation_z"] if pd.notna(r["conservation_z"]) else 0.0)
                banc_ok.append(r["BANC"])
            except ValueError:
                pass
    P = np.array(P) / 1000.0  # nm -> µm
    z = np.array(z)
    idx = {b: k for k, b in enumerate(banc_ok)}

    # conserved edges = BANC-induced edges among circuit nodes (all are consensus)
    def lf(fp, ns):
        d = pd.read_csv(fp, dtype=str); d.columns = ["s", "t"]
        d = d[d["s"].isin(ns) & d["t"].isin(ns)]
        return list(zip(d["s"], d["t"]))
    banc_edge = None
    for cand in ("banc_626_edge_list.csv", "banc_626_edge_list (2).csv"):
        if os.path.exists(dd + cand):
            banc_edge = dd + cand; break
    edges = lf(banc_edge, set(banc_ok)) if banc_edge else []

    norm = colors.Normalize(vmin=0, vmax=max(1.0, float(z.max())))
    cmap = cm.get_cmap("viridis")
    node_colors = [cmap(norm(v)) for v in z]

    fig = plt.figure(figsize=(8, 8)); fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111, projection="3d")

    def draw(angle):
        ax.clear()
        ax.scatter(P[:, 0], P[:, 1], P[:, 2], c=node_colors, s=40,
                   edgecolors="k", linewidths=0.3, depthshade=True)
        for u, v in edges:
            if u in idx and v in idx:
                a, b = idx[u], idx[v]
                ax.plot([P[a, 0], P[b, 0]], [P[a, 1], P[b, 1]],
                        [P[a, 2], P[b, 2]], color="#d4a017", lw=2.2)
        ax.set_title(f"Conserved circuit (N={len(P)}, {len(edges)} conserved edges)\n"
                     "colour = conservation z-score", fontsize=11)
        ax.set_axis_off()
        ax.view_init(elev=18, azim=angle)

    frames = np.linspace(0, 360, 60, endpoint=False)
    anim = FuncAnimation(fig, draw, frames=frames, interval=80)
    os.makedirs(figures_dir(), exist_ok=True)
    out = figures_dir() + "circuit_3d_conservation.gif"
    anim.save(out, writer=PillowWriter(fps=12), dpi=90)
    print(f"  Wrote {out}  ({len(P)} neurons, {len(edges)} conserved edges)")


if __name__ == "__main__":
    main()
