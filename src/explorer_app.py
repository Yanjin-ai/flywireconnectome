"""
Interactive explorer for the conserved circuit (Phase 1 / line C).
==================================================================

A single-page Streamlit app over the committed artifacts (no 300 MB data
download needed): browse the N=105 conserved circuit, filter by class, inspect
the per-neuron conservation track, view the conserved-edge subgraph, see the
enrichment / robustness numbers, and download a filtered CSV.

Run:  streamlit run src/explorer_app.py
"""
import json
import os

import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(name, default=None):
    p = os.path.join(REPO, "results", name)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return default


@st.cache_data
def load_data():
    enr = pd.read_csv(os.path.join(REPO, "network_enriched.csv"),
                      dtype={"BANC": str, "FAFB": str, "MANC": str})
    cons_p = os.path.join(REPO, "results", "neuron_conservation.csv")
    cons = (pd.read_csv(cons_p, dtype={"BANC": str})
            if os.path.exists(cons_p) else None)
    return enr, cons


CLASS_COLORS = {"descending": "#d62728", "ascending": "#1f77b4",
                "sensory_ascending": "#2ca02c", "sensory_descending": "#9467bd"}


def main():
    st.set_page_config(page_title="Conserved Connectome Circuit", layout="wide")
    st.title("Structural Invariance at the Sensorimotor Interface")
    st.caption("Conserved MCIS circuit across BANC × FAFB × MANC — interactive explorer")

    enr, cons = load_data()
    can = load("canonical_results.json", {})
    der = load("derived_stats.json", {})
    track = load("conservation_track.json", {})

    # headline metrics
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Circuit size N", can.get("N_reported", len(enr)))
    c2.metric("Conserved edges", can.get("n_conserved_edges", "—"))
    if track:
        c3.metric("Edge conservation vs degree-null",
                  f"{track['degree_null']['enrichment']:.1f}×",
                  f"Z={track['degree_null']['z']:.0f}σ")
    if der:
        c4.metric("Descending enrichment",
                  f"{der['enrichment']['descending']['fold']:.1f}×")

    st.sidebar.header("Filter")
    classes = sorted(enr["super_class"].dropna().unique())
    pick = st.sidebar.multiselect("Neuron class", classes, default=classes)
    view = enr[enr["super_class"].isin(pick)].copy()
    if cons is not None:
        view = view.merge(cons[["BANC", "conservation_z", "consensus_edges"]],
                          on="BANC", how="left")

    tab1, tab2, tab3 = st.tabs(["Conserved subgraph", "Conservation track",
                                "Composition & download"])

    with tab1:
        st.subheader("Conserved directed edges (identical in all three connectomes)")
        # rebuild conserved edges from the circuit if edge lists unavailable:
        st.info("The conserved-edge subgraph is rendered from network_enriched.csv "
                "cell-type labels. Conserved edges (gold) are the directed pairs "
                "present in BANC, FAFB and MANC.")
        # draw a simple graph of the filtered neurons + any known conserved edges
        edges = can.get("conserved_edge_celltypes", [])
        G = nx.DiGraph()
        for ct in view["cell_type"].dropna().unique():
            G.add_node(ct)
        for e in edges:
            G.add_edge(e[0], e[1])
        if G.number_of_nodes():
            fig, ax = plt.subplots(figsize=(9, 6))
            pos = nx.spring_layout(G, seed=1)
            cols = [CLASS_COLORS.get(
                view.loc[view["cell_type"] == n, "super_class"].iloc[0], "#aaa")
                if (view["cell_type"] == n).any() else "#aaa" for n in G.nodes()]
            nx.draw_networkx_nodes(G, pos, node_color=cols, node_size=120, ax=ax)
            nx.draw_networkx_edges(G, pos, edge_color="#d4a017", width=2,
                                   ax=ax, arrowsize=12)
            ax.axis("off")
            st.pyplot(fig)

    with tab2:
        st.subheader("Per-neuron conservation track (z vs degree-preserving null)")
        if cons is not None:
            top = view.sort_values("conservation_z", ascending=False).head(20)
            st.dataframe(top[["cell_type", "super_class", "consensus_edges",
                              "conservation_z"]], width="stretch")
            fig, ax = plt.subplots(figsize=(9, 4))
            z = view["conservation_z"].dropna().sort_values(ascending=False).values
            ax.plot(z, color="#1f77b4")
            ax.axhline(0, color="#999", lw=0.8)
            ax.set_xlabel("neuron rank"); ax.set_ylabel("conservation z")
            ax.spines[["top", "right"]].set_visible(False)
            st.pyplot(fig)
        else:
            st.warning("Run src/conservation_track.py to generate the track.")

    with tab3:
        st.subheader("Composition")
        comp = view["super_class"].value_counts()
        st.bar_chart(comp)
        nt = view["neurotransmitter_predicted"].value_counts()
        st.bar_chart(nt)
        st.download_button("Download filtered circuit CSV",
                           view.to_csv(index=False), "circuit_filtered.csv",
                           "text/csv")


if __name__ == "__main__":
    main()
