"""
One-page research summary for the challenge's Research Component (FAFB dataset).
================================================================================

Produces a single-page PDF/PNG covering the four required elements for ONE
dataset (FAFB):
  (1) network-graph visualisation of the conserved circuit,
  (2) a panel for the Codex 3D meshes (with the ready-to-paste neuron list +
      Codex link — Codex browsing is public, no token needed),
  (3) observations / biological hypothesis,
  (4) references.

Also writes results/codex_circuit_ids.txt (the 105 FAFB root IDs to paste into
codex.flywire.ai) and results/codex_links.md (instructions).

Run:  python src/make_research_summary.py
Output: research_summary_fafb.pdf / .png, results/codex_circuit_ids.txt,
        results/codex_links.md
"""
import json
import os
import sys

import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import repo_root, results_dir, figures_dir  # noqa

REPO = repo_root()
CODEX = "https://codex.flywire.ai/"


def load(name, default=None):
    p = os.path.join(results_dir(), name)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return default


def write_codex_ids(enr):
    os.makedirs(results_dir(), exist_ok=True)
    ids = enr["FAFB"].astype(str).tolist()
    with open(results_dir() + "codex_circuit_ids.txt", "w") as f:
        f.write("\n".join(ids) + "\n")
    with open(results_dir() + "codex_links.md", "w") as f:
        f.write(
            "# Codex 3D meshes — how to view the conserved circuit (no token)\n\n"
            f"The {len(ids)} FAFB root IDs of the conserved circuit are in "
            "`codex_circuit_ids.txt`.\n\n"
            f"1. Open Codex: <{CODEX}> (public, no login required to browse).\n"
            "2. Paste the root IDs (comma- or newline-separated) into the cell "
            "search / multi-cell view to load their 3D meshes.\n"
            "3. Rotate to the cervical-connective view and screenshot — that "
            "screenshot is the 'Codex 3D meshes' panel of the one-page summary.\n\n"
            "Alternatively open `results/neuroglancer_state.json` in "
            "<https://ngl.flywire.ai/> (Neuroglancer renders the same meshes, "
            "coloured by conservation z-score).\n\n"
            "Highlight the hub neurons that carry the 12 conserved edges: "
            "DNa15, DNg04, DNp58, DNp65, DNge076, DNge019, DNge020, DNp47, "
            "DNg79, DNp54, DNg02_g, AN09B033, AN06A027.\n"
        )
    return ids


def main():
    enr = pd.read_csv(os.path.join(REPO, "network_enriched.csv"),
                      dtype={"BANC": str, "FAFB": str, "MANC": str})
    can = load("canonical_results.json", {})
    der = load("derived_stats.json", {})
    track = load("conservation_track.json", {})
    ids = write_codex_ids(enr)

    n = len(enr)
    desc = der.get("enrichment", {}).get("descending", {})
    asc = der.get("enrichment", {}).get("ascending", {})
    nz = track.get("degree_null", {})

    fig = plt.figure(figsize=(8.27, 11.69)); fig.patch.set_facecolor("white")
    fig.text(0.5, 0.975, "A Conserved Sensorimotor Circuit across Three "
             "Drosophila Connectomes", ha="center", va="top", fontsize=13,
             fontweight="bold")
    fig.text(0.5, 0.952, "Research summary · dataset: FAFB (♀ adult brain, "
             "Dorkenwald et al. 2024) · Yanjin Li", ha="center", fontsize=9,
             color="#444")

    # (1) network graph
    axg = fig.add_axes([0.06, 0.66, 0.44, 0.25]); axg.axis("off")
    axg.set_title("(1) Conserved-circuit network graph", fontsize=10,
                  fontweight="bold", loc="left")
    edges = can.get("conserved_edge_celltypes", [])
    G = nx.DiGraph(); G.add_edges_from([tuple(e) for e in edges])
    if G.number_of_nodes():
        pos = nx.spring_layout(G, seed=2, k=0.9)
        nx.draw_networkx_nodes(G, pos, node_color="#d62728", node_size=260, ax=axg)
        nx.draw_networkx_edges(G, pos, edge_color="#d4a017", width=2,
                               arrowsize=12, ax=axg)
        nx.draw_networkx_labels(G, pos, font_size=6, ax=axg)
    axg.text(0, -0.05, f"{n} neurons · {can.get('n_conserved_edges','?')} "
             "directed edges identical in BANC×FAFB×MANC (gold).",
             transform=axg.transAxes, fontsize=7, color="#444")

    # (2) Codex 3D meshes panel (instructions + anatomical render stand-in)
    axm = fig.add_axes([0.52, 0.66, 0.44, 0.25]); axm.axis("off")
    axm.set_title("(2) Codex 3D meshes", fontsize=10, fontweight="bold", loc="left")
    f6 = os.path.join(figures_dir(), "figure6_spatial.png")
    if os.path.exists(f6):
        axm.imshow(mpimg.imread(f6))
    axm.text(0, -0.06, f"Paste results/codex_circuit_ids.txt ({n} FAFB root IDs) "
             f"into {CODEX} → 3D meshes (public, no token). Anatomical "
             "distribution shown; meshes concentrate on the cervical connective.",
             transform=axm.transAxes, fontsize=7, color="#444", wrap=True)

    # (3) observations / hypothesis
    obs = (
        "(3) Observations & biological hypothesis\n\n"
        f"• The conserved set is almost entirely sensorimotor: descending "
        f"{desc.get('fold','?')}× and ascending {asc.get('fold','?')}× enriched "
        "vs the FAFB whole-brain background (DN+AN ≈ 93% of the circuit), giving "
        "direct structural support for the sensorimotor bottleneck / effectome "
        "(Pospisil et al. 2024).\n"
        f"• Specific wiring is conserved far beyond degree sequence: "
        f"{track.get('observed_consensus_edges','?')} all-three consensus edges "
        f"vs {nz.get('mean',0):.0f} expected under a degree-preserving null "
        f"({nz.get('enrichment',0):.1f}×, Z={nz.get('z',0):.0f}σ) — the agreement "
        "reflects real connectivity identity, not matched degree.\n"
        f"• {der.get('sexual_conservation',{}).get('pct_isomorphic','?')}% of the "
        "circuit is conserved across sexes (♀ FAFB/BANC vs ♂ MANC); the few "
        "dimorphic neurons project to abdominal VNC / lateral brain — the classes "
        "expected to diverge for reproductive behaviour.\n"
        "• The conserved edges form reciprocal DN↔DN pairs (DNa15↔DNg04, "
        "DNp58↔DNp65) with mixed ACh/GABA chemistry — a feedforward-inhibition "
        "motif (Milo et al. 2002) plausibly filtering descending motor commands.\n"
        "• Hypothesis: these edges are a developmentally canalised brain↔cord "
        "communication backbone; silencing the hub neurons should impair walking, "
        "flight and posture simultaneously (testable via optogenetics)."
    )
    fig.text(0.06, 0.61, obs, fontsize=8.4, va="top", ha="left", wrap=True)

    refs = (
        "(4) References\n"
        "1. Dorkenwald et al. (2024) Nature 634:123 — FAFB connectome.\n"
        "2. Schlegel et al. (2024) Nature 634:139 — multi-connectome cell typing / NBLAST.\n"
        "3. Bates et al. (2025) bioRxiv — BANC; cross-dataset correspondence metadata.\n"
        "4. Pospisil et al. (2024) Nature 634:234 — sensorimotor bottleneck / effectome.\n"
        "5. Witvliet et al. (2021) Nature 596:257 — connectome stereotypy.\n"
        "6. Milo et al. (2002) Science 298:824 — network motifs (feedforward inhibition)."
    )
    fig.text(0.06, 0.20, refs, fontsize=8, va="top", ha="left")
    fig.text(0.06, 0.03, "Code & reproducible results: "
             "github.com/Yanjin-ai/flywireconnectome · live demo: "
             "yanjin-ai-flywireconnectome-srcexplorer-app-pmdboy.streamlit.app",
             fontsize=7, color="#666")

    pdf = os.path.join(REPO, "research_summary_fafb.pdf")
    with PdfPages(pdf) as pp:
        pp.savefig(fig, facecolor="white")
    fig.savefig(os.path.join(REPO, "research_summary_fafb.png"), dpi=130,
                facecolor="white")
    plt.close(fig)
    print(f"  Wrote {pdf} (1 page) + research_summary_fafb.png")
    print(f"  Wrote results/codex_circuit_ids.txt ({len(ids)} FAFB IDs) + "
          "results/codex_links.md")


if __name__ == "__main__":
    main()
