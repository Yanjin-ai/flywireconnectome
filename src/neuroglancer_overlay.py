"""
Neuroglancer overlay for the conserved circuit (Phase 1 / line C).
==================================================================

Generates a FlyWire-compatible Neuroglancer state (JSON) that loads the
conserved circuit's FAFB neurons, coloured either by neuron class or by the
per-neuron conservation z-score (the "conservation track"). This is the
ecosystem-native deliverable: the state can be opened in the FlyWire
Neuroglancer / shared as a short link.

Dependency-free: the state is plain JSON. If `fafbseg` is installed and you are
authenticated, pass --shorten to also produce a shortened FlyWire share URL via
fafbseg.flywire.encode_url; otherwise the full JSON state is written and the
manual sharing instructions are printed.

Run:  python src/neuroglancer_overlay.py [--color class|conservation] [--shorten]
Output: results/neuroglancer_state.json
"""
import argparse
import json
import os
import sys

import pandas as pd
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import repo_root, results_dir  # noqa

# Public FlyWire FAFB flat segmentation (materialisation 783). Adjust if you
# target a different version; the circuit's FAFB column uses v783 root IDs.
FLYWIRE_SEG_SOURCE = "precomputed://gs://flywire_v141_m783"
FLYWIRE_NG = "https://ngl.flywire.ai/"

CLASS_COLORS = {
    "descending": "#d62728", "ascending": "#1f77b4",
    "sensory_ascending": "#2ca02c", "sensory_descending": "#9467bd",
}


def conservation_colors(zs):
    norm = mcolors.Normalize(vmin=0, vmax=max(1.0, float(zs.max())))
    cmap = cm.get_cmap("viridis")
    return [mcolors.to_hex(cmap(norm(z))) for z in zs]


def build_state(df, color_by):
    segs = df["FAFB"].astype(str).tolist()
    if color_by == "conservation" and "conservation_z" in df:
        cols = conservation_colors(df["conservation_z"].fillna(0).values)
    else:
        cols = [CLASS_COLORS.get(c, "#aaaaaa") for c in df["super_class"]]
    seg_colors = {s: c for s, c in zip(segs, cols)}
    state = {
        "dimensions": {"x": [4e-9, "m"], "y": [4e-9, "m"], "z": [4e-8, "m"]},
        "position": [120000, 50000, 2000],
        "crossSectionScale": 2.0,
        "projectionScale": 60000,
        "layers": [
            {
                "type": "segmentation",
                "source": FLYWIRE_SEG_SOURCE,
                "tab": "segments",
                "segments": segs,
                "segmentColors": seg_colors,
                "name": f"conserved-circuit ({color_by})",
            }
        ],
        "selectedLayer": {"layer": f"conserved-circuit ({color_by})",
                          "visible": True},
        "layout": "3d",
    }
    return state


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--color", choices=["class", "conservation"], default="class")
    ap.add_argument("--shorten", action="store_true")
    args = ap.parse_args()

    repo = repo_root()
    df = pd.read_csv(os.path.join(repo, "network_enriched.csv"),
                     dtype={"BANC": str, "FAFB": str})
    cons_path = os.path.join(results_dir(), "neuron_conservation.csv")
    if os.path.exists(cons_path):
        cons = pd.read_csv(cons_path, dtype={"BANC": str})
        df = df.merge(cons[["BANC", "conservation_z"]], on="BANC", how="left")

    state = build_state(df, args.color)
    os.makedirs(results_dir(), exist_ok=True)
    out = os.path.join(results_dir(), "neuroglancer_state.json")
    with open(out, "w") as f:
        json.dump(state, f, indent=2)
    print(f"  Wrote {out}  ({len(df)} segments, coloured by {args.color})")
    print(f"  Segmentation source: {FLYWIRE_SEG_SOURCE}")

    if args.shorten:
        try:
            from fafbseg.flywire import encode_url
            url = encode_url(segments=df["FAFB"].astype(int).tolist(),
                             seg_colors={int(k): v
                                         for k, v in state["layers"][0]
                                         ["segmentColors"].items()},
                             open_browser=False)
            print(f"  Shortened FlyWire URL: {url}")
        except Exception as e:
            print(f"  (--shorten unavailable: {e})")
            print(f"  Open {FLYWIRE_NG} and load the JSON state manually, or use "
                  "fafbseg.flywire.encode_url when authenticated.")
    else:
        print(f"  To share: open {FLYWIRE_NG} → paste the JSON state, or run "
              "with --shorten (needs fafbseg + FlyWire auth).")


if __name__ == "__main__":
    main()
