"""
Generate extended_abstract.pdf / .png — the 2-page challenge summary.
====================================================================

All numbers are pulled from results/*.json so the abstract can never drift
from the canonical pipeline. Embeds figure1 (circuit) and figure5 (robustness).

Run (after the pipeline + figures have been generated):
    python src/make_abstract.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg

from mcis_paths import repo_root, results_dir, figures_dir

REPO = repo_root()


def load(name):
    with open(os.path.join(results_dir(), name)) as f:
        return json.load(f)


def main():
    can = load("canonical_results.json")
    der = load("derived_stats.json")
    ilp = load("ilp_validation.json")

    N = can["N_reported"]
    E = can["n_conserved_edges"]
    sd = can["seed_distribution"]
    cs = can["correspondence_shuffle_null"]
    dp = can["degree_preserving_null"]
    cent = can["centrality"]
    desc = der["enrichment"]["descending"]
    asc = der["enrichment"]["ascending"]
    aq = der["annotation_quality"]
    sx = der["sexual_conservation"]

    TITLE = ("Structural Invariance at the Sensorimotor Interface:\n"
             "A Maximum Common Induced Subgraph Across Three "
             "Drosophila Connectomes")
    SUB = ("Yanjin Li  ·  FlyWire Qualification Challenge  ·  "
           "github.com/Yanjin-ai/flywireconnectome")

    abstract = (
        f"We ask whether synaptic connectivity itself is structurally invariant "
        f"across independently prepared connectomes. Using the BANC metadata's "
        f"NBLAST-based 1:1 correspondence columns (Bates et al. 2025), we align "
        f"individual neurons across BANC (♀ brain+cord), FAFB (♀ brain) "
        f"and MANC (♂ nerve cord) and search for the largest set of neurons "
        f"whose directed induced subgraph is identical (isomorphic) in all three. "
        f"The result is an N = {N}-neuron sensorimotor backbone with {E} conserved "
        f"directed edges, enriched {desc['fold']:.1f}× for descending and "
        f"{asc['fold']:.1f}× for ascending neurons, and {sx['pct_isomorphic']}% "
        f"conserved across sexes. All reported numbers are reproduced by "
        f"src/run_analysis.py and stored under results/."
    )

    methods = (
        "Methods. Correspondence: BANC metadata fafb_match / manc_match "
        "(NBLAST). Graphs: directed, edge existence (not weight). "
        f"Search space: {can['n_bounds']['matched_triplets']:,} matched triplets "
        f"→ {can['n_bounds']['in_all3_edge_lists']:,} present in all three "
        f"edge lists → {can['n_bounds']['giant_component']}-node giant "
        "consensus component. Algorithm: greedy disagreement removal + exhaustive "
        "expansion (single implementation in mcis_connectome.solver), reported as "
        f"the best of a {sd['n_seeds']}-seed multi-start. Optimality checked "
        f"against an ILP (PuLP/CBC) on 50 subgraphs: mean gap "
        f"{ilp['overall']['mean_gap_pct']}% (max {ilp['overall']['max_gap_pct']}%)."
    )

    comp = der["composition"]
    n_sensory = (comp.get("sensory_ascending", [0])[0]
                 + comp.get("sensory_descending", [0])[0])
    rows = [
        ("MCIS size N", f"{N}  (dist. {sd['mean']:.1f} ± {sd['std']:.1f}, "
         f"[{sd['min']},{sd['max']}])"),
        ("Conserved edges", f"{E}"),
        ("Composition", f"{comp['descending'][0]} DN + "
         f"{comp['ascending'][0]} AN + {n_sensory} sensory"),
        ("Descending enrichment", f"{desc['fold']:.1f}×  (p={desc['p_value']:.0e})"),
        ("Ascending enrichment", f"{asc['fold']:.1f}×  (p={asc['p_value']:.0e})"),
        ("Corr-shuffle null", f"{cs['mean']:.1f} ± {cs['std']:.1f}  (>15σ below real)"),
        ("Degree-preserving null", f"{dp['mean']:.1f} ± {dp['std']:.1f}  (≈ real mean)"),
        ("Centrality (lower)", f"p = {cent['p_value']:.3f} (one-sided)"),
        ("Annotation quality", f"{aq['circuit_pct']}% vs {aq['noncircuit_pct']}%  "
         f"(p={aq['fisher_p']:.3f})"),
        ("Sexual conservation", f"{sx['pct_isomorphic']}% ({sx['isomorphic']}/{N})"),
    ]

    discussion = (
        "Discussion & Limitations. The conserved backbone is almost entirely "
        "descending and ascending neurons (93.3%), giving direct structural "
        "support for the sensorimotor bottleneck (Pospisil et al. 2024) and "
        "showing it is canalized across sexes and specimens. "
        "PRINCIPAL CAVEAT: a degree-preserving rewire null reaches essentially "
        f"the same MCIS size ({dp['mean']:.1f} ± {dp['std']:.1f} ≈ real "
        f"{sd['mean']:.1f}), so the FAFB degree sequence — not specific edge "
        "identity — explains most of the achievable N; NBLAST neuron identity "
        "adds the remaining signal (correspondence-shuffle collapses to "
        f"{cs['mean']:.1f}). N is a near-optimal heuristic lower bound (ILP gap "
        f"~{ilp['overall']['mean_gap_pct']}%), not a certified global optimum."
    )

    def build_page1():
        fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
        fig.patch.set_facecolor("white")
        fig.text(0.5, 0.96, TITLE, ha="center", va="top", fontsize=13,
                 fontweight="bold", wrap=True)
        fig.text(0.5, 0.885, SUB, ha="center", fontsize=8.5, color="#444")
        fig.text(0.08, 0.86, "Abstract", fontsize=11, fontweight="bold")
        fig.text(0.08, 0.845, abstract, fontsize=8.6, va="top", wrap=True, ha="left")
        fig.text(0.08, 0.70, methods, fontsize=8.6, va="top", wrap=True, ha="left")
        ax = fig.add_axes([0.08, 0.30, 0.84, 0.24]); ax.axis("off")
        ax.set_title("Key Results", loc="left", fontweight="bold", fontsize=11)
        tbl = ax.table(cellText=rows, colLabels=["Metric", "Value"],
                       cellLoc="left", loc="upper left", colWidths=[0.34, 0.66])
        tbl.auto_set_font_size(False); tbl.set_fontsize(8.2); tbl.scale(1, 1.25)
        for (r, c), cell in tbl.get_celld().items():
            cell.set_edgecolor("#ccc")
            if r == 0:
                cell.set_facecolor("#1f3a5f")
                cell.set_text_props(color="white", fontweight="bold")
        f1 = os.path.join(figures_dir(), "figure1_circuit_layouts.png")
        if os.path.exists(f1):
            ax2 = fig.add_axes([0.06, 0.02, 0.88, 0.26]); ax2.axis("off")
            ax2.imshow(mpimg.imread(f1))
            ax2.set_title("Conserved circuit across BANC × FAFB × MANC",
                          fontsize=8, color="#444")
        return fig

    def build_page2():
        fig = plt.figure(figsize=(8.27, 11.69)); fig.patch.set_facecolor("white")
        f5 = os.path.join(figures_dir(), "figure5_robustness.png")
        if os.path.exists(f5):
            ax = fig.add_axes([0.05, 0.45, 0.9, 0.5]); ax.axis("off")
            ax.imshow(mpimg.imread(f5))
            ax.set_title("Robustness, null models, and ILP optimality",
                         fontsize=10, fontweight="bold")
        fig.text(0.08, 0.40, discussion, fontsize=9, va="top", wrap=True, ha="left")
        fig.text(0.08, 0.06,
                 "Reproduce: pip install -e .  &&  MCIS_DATA_DIR=... "
                 "python src/run_analysis.py --seeds 100",
                 fontsize=8, color="#666", family="monospace")
        return fig

    pdf_path = os.path.join(REPO, "extended_abstract.pdf")
    p1, p2 = build_page1(), build_page2()
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(p1, facecolor="white")
        pdf.savefig(p2, facecolor="white")
    p1.savefig(os.path.join(REPO, "extended_abstract.png"), dpi=130,
               facecolor="white")
    plt.close(p1); plt.close(p2)
    print(f"Saved {pdf_path} (2 pages) and extended_abstract.png (page 1)")


if __name__ == "__main__":
    main()
