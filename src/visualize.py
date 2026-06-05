"""
Circuit Visualization — 75-neuron conserved sensorimotor circuit
BANC × FAFB × MANC
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = '/Volumes/SANDISK ELE/flywire研究/'
OUT_DIR  = '/Volumes/SANDISK ELE/flywire研究/repo/figures/'

# ── load data ────────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_DIR + 'network_enriched.csv')
net = pd.read_csv(DATA_DIR + 'network.csv')

# Rebuild the verified induced subgraph
def load_filtered(fp, nodes):
    d = pd.read_csv(fp, dtype=str); d.columns = ['s','t']
    d = d[d['s'].isin(nodes) & d['t'].isin(nodes)]
    return nx.from_pandas_edgelist(d, 's', 't', create_using=nx.DiGraph())

banc_nodes = set(net['BANC'].astype(str))
fafb_nodes = set(net['FAFB'].astype(str))
manc_nodes = set(net['MANC'].astype(str))

Gb = load_filtered(DATA_DIR + 'banc_626_edge_list (2).csv', banc_nodes)
Gf = load_filtered(DATA_DIR + 'fafb_783_edge_list.csv', fafb_nodes)
Gm = load_filtered(DATA_DIR + 'manc_1.2.1_edge_list.csv', manc_nodes)

# Map BANC→indices
banc_list = list(net['BANC'].astype(str))
bi = {b: i for i, b in enumerate(banc_list)}
fi = {f: i for i, f in enumerate(net['FAFB'].astype(str))}
mi = {m: i for i, m in enumerate(net['MANC'].astype(str))}

be = {(bi[u], bi[v]) for u,v in Gb.edges() if u in bi and v in bi}
fe = {(fi[u], fi[v]) for u,v in Gf.edges() if u in fi and v in fi}
me = {(mi[u], mi[v]) for u,v in Gm.edges() if u in mi and v in mi}
circuit_edges = be & fe & me   # 6 verified edges

# Build circuit graph (BANC IDs as labels)
CG = nx.DiGraph()
for idx, row in df.iterrows():
    CG.add_node(idx,
                cell_type  = str(row.get('cell_type', '?')),
                super_class= str(row.get('super_class', 'unknown')),
                nt         = str(row.get('neurotransmitter_predicted', 'unknown')),
                cns_net    = str(row.get('cns_network', 'unknown')),
                dimorphic  = str(row.get('sexually_dimorphic', 'unknown')))
for u_idx, v_idx in circuit_edges:
    CG.add_edge(u_idx, v_idx)

# ── colour maps ───────────────────────────────────────────────────────────────
SC_COLOR = {
    'descending':        '#E63946',
    'ascending':         '#2196F3',
    'sensory_ascending': '#4CAF50',
    'sensory_descending':'#FF9800',
    'unknown':           '#9E9E9E',
}
NT_COLOR = {
    'acetylcholine': '#4FC3F7',
    'gaba':          '#EF9A9A',
    'glutamate':     '#A5D6A7',
    'serotonin':     '#CE93D8',
    'dopamine':      '#FFCC02',
    'unknown':       '#BDBDBD',
}
CNS_COLOR = {
    'leg VNC':              '#1565C0',
    'dorsal VNC':           '#283593',
    'flange median bundle': '#0277BD',
    'abdominal VNC':        '#00695C',
    'lateral brain':        '#558B2F',
    'posterior brain':      '#827717',
    'unknown':              '#616161',
}

def sc_color(n):  return SC_COLOR.get(CG.nodes[n]['super_class'], '#9E9E9E')
def nt_color(n):  return NT_COLOR.get(CG.nodes[n]['nt'], '#BDBDBD')
def cns_color(n): return CNS_COLOR.get(CG.nodes[n]['cns_net'], '#616161')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Main circuit diagram
# ══════════════════════════════════════════════════════════════════════════════
fig1, axes = plt.subplots(1, 3, figsize=(21, 8))
fig1.patch.set_facecolor('#0D1117')

layouts = {
    'spring':  nx.spring_layout(CG, seed=42, k=1.8),
    'kamada':  nx.kamada_kawai_layout(CG),
    'circular':nx.circular_layout(CG),
}
layout_names = ['spring', 'kamada', 'circular']
subtitles = ['Spring Layout', 'Kamada–Kawai Layout', 'Circular Layout']

for ax, lname, subtitle in zip(axes, layout_names, subtitles):
    ax.set_facecolor('#0D1117')
    pos = layouts[lname]
    node_colors = [sc_color(n) for n in CG.nodes()]
    # Highlight nodes with edges
    node_sizes  = [200 if CG.degree(n) > 0 else 60 for n in CG.nodes()]
    node_border = ['white' if CG.degree(n) > 0 else 'none' for n in CG.nodes()]

    nx.draw_networkx_nodes(CG, pos, ax=ax,
                           node_color=node_colors,
                           node_size=node_sizes,
                           edgecolors=node_border,
                           linewidths=0.8, alpha=0.9)
    nx.draw_networkx_edges(CG, pos, ax=ax,
                           edge_color='#FFD700',
                           arrows=True,
                           arrowsize=18,
                           width=2.5,
                           connectionstyle='arc3,rad=0.15',
                           alpha=0.95)
    # Label only the 6 hub nodes
    hub_nodes = [n for n in CG.nodes() if CG.degree(n) > 0]
    labels = {n: CG.nodes[n]['cell_type'][:10] for n in hub_nodes}
    nx.draw_networkx_labels(CG, pos, labels, ax=ax,
                            font_size=6.5, font_color='white', font_weight='bold')
    ax.set_title(subtitle, color='white', fontsize=13, pad=10)
    ax.axis('off')

# Legend
legend_elements = [
    mpatches.Patch(facecolor=v, label=k.replace('_',' ').title())
    for k,v in SC_COLOR.items() if k != 'unknown'
]
axes[0].legend(handles=legend_elements, loc='lower left',
               framealpha=0.25, labelcolor='white',
               facecolor='#1a1a2e', fontsize=8, title='Neuron class',
               title_fontsize=9)
legend_elements[0].set_label('Descending (DN)')
legend_elements[1].set_label('Ascending (AN)')
legend_elements[2].set_label('Sensory-ascending')
legend_elements[3].set_label('Sensory-descending')
axes[0].legend(handles=legend_elements, loc='lower left',
               framealpha=0.25, labelcolor='white',
               facecolor='#1a1a2e', fontsize=8, title='Neuron class',
               title_fontsize=9)

fig1.suptitle(
    'Conserved Sensorimotor Circuit  ·  N = 75 neurons  ·  BANC × FAFB × MANC',
    color='white', fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
fig1.savefig(OUT_DIR + 'figure1_circuit_layouts.png',
             dpi=180, bbox_inches='tight', facecolor='#0D1117')
print('Saved figure1_circuit_layouts.png')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Composition + connectivity breakdown
# ══════════════════════════════════════════════════════════════════════════════
fig2 = plt.figure(figsize=(20, 12))
fig2.patch.set_facecolor('#0D1117')
gs = gridspec.GridSpec(2, 3, figure=fig2, hspace=0.45, wspace=0.35)

ax_sc   = fig2.add_subplot(gs[0, 0])
ax_nt   = fig2.add_subplot(gs[0, 1])
ax_cns  = fig2.add_subplot(gs[0, 2])
ax_deg  = fig2.add_subplot(gs[1, 0])
ax_mat  = fig2.add_subplot(gs[1, 1:])

def dark_bar(ax, labels, values, colors, title, ylabel='Count'):
    ax.set_facecolor('#161B22')
    bars = ax.bar(labels, values, color=colors, width=0.6, edgecolor='none', alpha=0.9)
    ax.set_title(title, color='white', fontsize=11, pad=8)
    ax.set_ylabel(ylabel, color='#aaa', fontsize=9)
    ax.tick_params(colors='#aaa', labelsize=8)
    ax.spines[:].set_color('#333')
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                str(v), ha='center', va='bottom', color='white', fontsize=8, fontweight='bold')
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right', color='#ccc')

# Superclass bar
sc_vc = df['super_class'].value_counts()
dark_bar(ax_sc,
         [s.replace('_','\n') for s in sc_vc.index],
         sc_vc.values,
         [SC_COLOR.get(s,'#9E9E9E') for s in sc_vc.index],
         'A  Neuron Class Composition')

# Neurotransmitter bar
nt_vc = df['neurotransmitter_predicted'].value_counts()
dark_bar(ax_nt,
         nt_vc.index.tolist(),
         nt_vc.values,
         [NT_COLOR.get(n,'#BDBDBD') for n in nt_vc.index],
         'B  Neurotransmitter Profile')

# CNS target network bar
cns_vc = df['cns_network'].value_counts().head(7)
dark_bar(ax_cns,
         [c.replace(' ','\n') for c in cns_vc.index],
         cns_vc.values,
         [CNS_COLOR.get(c,'#616161') for c in cns_vc.index],
         'C  Motor Target Region')

# Degree distribution
ax_deg.set_facecolor('#161B22')
degrees = [CG.degree(n) for n in CG.nodes()]
deg_vals, deg_counts = np.unique(degrees, return_counts=True)
ax_deg.bar(deg_vals, deg_counts, color='#E63946', alpha=0.85, edgecolor='none')
ax_deg.set_title('D  Node Degree Distribution', color='white', fontsize=11, pad=8)
ax_deg.set_xlabel('Degree', color='#aaa', fontsize=9)
ax_deg.set_ylabel('Count', color='#aaa', fontsize=9)
ax_deg.tick_params(colors='#aaa')
ax_deg.spines[:].set_color('#333')

# Connectivity comparison matrix across datasets
ax_mat.set_facecolor('#161B22')
datasets = ['BANC', 'FAFB', 'MANC']
edge_counts = [len(be), len(fe), len(me)]
consensus = len(circuit_edges)
x = np.arange(3)
bars = ax_mat.bar(x, edge_counts, color=['#4FC3F7','#A5D6A7','#EF9A9A'],
                  width=0.5, alpha=0.85, label='Dataset edges')
ax_mat.axhline(consensus, color='#FFD700', lw=2.5, ls='--', label=f'Consensus ({consensus} edges)')
ax_mat.set_xticks(x); ax_mat.set_xticklabels(datasets, color='#ccc', fontsize=11)
ax_mat.set_title('E  Edge Count Comparison — Same N=75 Neurons Across Datasets',
                 color='white', fontsize=11, pad=8)
ax_mat.set_ylabel('Synapse Connections', color='#aaa', fontsize=9)
ax_mat.tick_params(colors='#aaa')
ax_mat.spines[:].set_color('#333')
ax_mat.legend(facecolor='#1a1a2e', labelcolor='white', fontsize=9)
for bar, v in zip(bars, edge_counts):
    ax_mat.text(bar.get_x()+bar.get_width()/2, bar.get_height()+200,
                f'{v:,}', ha='center', color='white', fontsize=10, fontweight='bold')

fig2.suptitle(
    'Circuit Composition & Cross-Connectome Consistency  ·  75-neuron Shared Circuit',
    color='white', fontsize=14, fontweight='bold', y=1.01)
fig2.savefig(OUT_DIR + 'figure2_composition.png',
             dpi=180, bbox_inches='tight', facecolor='#0D1117')
print('Saved figure2_composition.png')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Hub neurons + 6 conserved edges highlighted
# ══════════════════════════════════════════════════════════════════════════════
fig3, ax = plt.subplots(figsize=(14, 10))
fig3.patch.set_facecolor('#0D1117')
ax.set_facecolor('#0D1117')

pos = nx.spring_layout(CG, seed=42, k=2.5)

# Draw non-hub nodes (small, dim)
non_hub = [n for n in CG.nodes() if CG.degree(n) == 0]
hub     = [n for n in CG.nodes() if CG.degree(n)  > 0]
nx.draw_networkx_nodes(CG, pos, nodelist=non_hub, ax=ax,
                       node_color=[sc_color(n) for n in non_hub],
                       node_size=35, alpha=0.25)
nx.draw_networkx_nodes(CG, pos, nodelist=hub, ax=ax,
                       node_color=[sc_color(n) for n in hub],
                       node_size=350, edgecolors='white', linewidths=1.5, alpha=0.95)

# Draw edges with glow effect (multiple passes)
for width, alpha in [(6, 0.10), (4, 0.20), (2.5, 0.95)]:
    nx.draw_networkx_edges(CG, pos, ax=ax,
                           edge_color='#FFD700', arrows=True,
                           arrowsize=20, width=width,
                           connectionstyle='arc3,rad=0.2', alpha=alpha)

# Label hub nodes
labels = {n: CG.nodes[n]['cell_type'] for n in hub}
nx.draw_networkx_labels(CG, pos, labels, ax=ax,
                        font_size=9, font_color='white', font_weight='bold')

# Annotate edge with NT info
for (u, v) in CG.edges():
    ux, uy = pos[u]
    vx, vy = pos[v]
    mx, my = (ux+vx)/2, (uy+vy)/2
    nt_u = CG.nodes[u]['nt']
    ax.annotate(f'{nt_u[:3]}', xy=(mx, my), fontsize=7.5, color='#FFD700',
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.2', fc='#0D1117', alpha=0.6))

# Legend
sc_patches = [mpatches.Patch(color=v, label=k.replace('_',' ').title())
              for k, v in SC_COLOR.items() if k != 'unknown']
ax.legend(handles=sc_patches, loc='lower left',
          facecolor='#161B22', labelcolor='white', fontsize=9,
          title='Neuron class', title_fontsize=10, framealpha=0.8)

ax.set_title(
    '6 Conserved Synaptic Connections  ·  Verified across BANC × FAFB × MANC\n'
    'Gold = conserved edge  ·  Large nodes = circuit hubs  ·  Small nodes = isolated members',
    color='white', fontsize=12, pad=12)
ax.axis('off')
plt.tight_layout()
fig3.savefig(OUT_DIR + 'figure3_hub_circuit.png',
             dpi=180, bbox_inches='tight', facecolor='#0D1117')
print('Saved figure3_hub_circuit.png')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Sexual dimorphism + neurotransmitter scatter
# ══════════════════════════════════════════════════════════════════════════════
fig4, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig4.patch.set_facecolor('#0D1117')

# Dimorphism pie
ax1.set_facecolor('#0D1117')
dim_vc = df['sexually_dimorphic'].value_counts()
colors_pie = ['#4FC3F7' if 'iso' in str(k).lower() else '#EF9A9A' for k in dim_vc.index]
wedges, texts, autotexts = ax1.pie(
    dim_vc.values, labels=dim_vc.index,
    colors=colors_pie, autopct='%1.0f%%',
    startangle=90, pctdistance=0.75,
    wedgeprops=dict(linewidth=2, edgecolor='#0D1117'))
for t in texts:     t.set_color('white')
for t in autotexts: t.set_color('white'); t.set_fontsize(11); t.set_fontweight('bold')
ax1.set_title('F  Sexual Dimorphism\n(isomorphic = conserved across sexes)',
              color='white', fontsize=11, pad=10)

# NT × superclass scatter (stacked bar)
ax2.set_facecolor('#161B22')
sc_nt = df.groupby(['super_class','neurotransmitter_predicted']).size().unstack(fill_value=0)
sc_order = ['descending','ascending','sensory_ascending','sensory_descending']
sc_nt = sc_nt.reindex([s for s in sc_order if s in sc_nt.index])
bottom = np.zeros(len(sc_nt))
for nt_name in sc_nt.columns:
    ax2.bar(range(len(sc_nt)), sc_nt[nt_name], bottom=bottom,
            color=NT_COLOR.get(nt_name,'#BDBDBD'), label=nt_name, alpha=0.88)
    bottom += sc_nt[nt_name].values
ax2.set_xticks(range(len(sc_nt)))
ax2.set_xticklabels([s.replace('_','\n') for s in sc_nt.index], color='#ccc', fontsize=9)
ax2.set_title('G  Neurotransmitter × Neuron Class', color='white', fontsize=11, pad=8)
ax2.set_ylabel('Neuron count', color='#aaa', fontsize=9)
ax2.tick_params(colors='#aaa')
ax2.spines[:].set_color('#333')
ax2.legend(facecolor='#1a1a2e', labelcolor='white', fontsize=9, loc='upper right')

fig4.suptitle('Sexual Conservation & Neurochemical Profile', color='white', fontsize=13, y=1.01)
plt.tight_layout()
fig4.savefig(OUT_DIR + 'figure4_dimorphism_nt.png',
             dpi=180, bbox_inches='tight', facecolor='#0D1117')
print('Saved figure4_dimorphism_nt.png')

print('\nAll figures saved to', OUT_DIR)
