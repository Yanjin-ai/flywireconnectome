"""
Regenerate figure7_nblast_confidence.png with a clean white background.

Data from science.md §4.5 (NBLAST confidence curve) and §2.1 (triplet pool).
Run from repo root:  python src/regenerate_figure7.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT = 'figures/figure7_nblast_confidence.png'

# ── Data ────────────────────────────────────────────────────────────────────
# Panel A: NBLAST confidence tier → MCIS size (science.md §4.5)
tiers    = [10, 20, 30, 50, 75, 100]          # top-k %
pools    = [279, 559, 839, 1399, 2098, 2798]  # triplet pool at each tier
mcis_n   = [10,  28,  30,  59,  88,  104]     # MCIS N at each tier

# Panel B: NBLAST agreement distribution across 2,798 triplets
agree_labels = ['Neither agrees\n(expert override)', 'One agrees\n(partial)', 'Both agree\n(high confidence)']
agree_counts = [1024, 1380, 394]   # total = 2798
agree_pcts   = [f'{v/2798:.0%}' for v in agree_counts]
agree_colors = ['#d62728', '#ff7f0e', '#2ca02c']   # red / orange / green

# ── Figure layout ────────────────────────────────────────────────────────────
fig, (ax_a, ax_b) = plt.subplots(
    1, 2,
    figsize=(13, 5.2),
    gridspec_kw={'width_ratios': [1.4, 1]},
)
fig.patch.set_facecolor('white')
fig.suptitle(
    'NBLAST Confidence Analysis — Matching Quality vs Circuit Conservation',
    fontsize=13, fontweight='bold', y=1.01, color='#1a1a1a',
)

# ── Panel A: confidence curve ────────────────────────────────────────────────
ax_a.set_facecolor('#f8f9fa')
ax_a.fill_between(tiers, mcis_n, alpha=0.18, color='#1f77b4', zorder=1)
ax_a.plot(tiers, mcis_n,
          color='#1f77b4', linewidth=2.4, marker='o',
          markersize=8, markerfacecolor='white', markeredgewidth=2, zorder=3)

for x, y, pool in zip(tiers, mcis_n, pools):
    ax_a.annotate(
        f'N={y}\n({pool:,} pool)',
        xy=(x, y),
        xytext=(0, 14),
        textcoords='offset points',
        ha='center', va='bottom',
        fontsize=8.5, color='#1f77b4',
        arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.8),
    )

ax_a.set_xlabel('Top-k% of triplets (ranked by NBLAST agreement)', fontsize=10)
ax_a.set_ylabel('MCIS size N', fontsize=10)
ax_a.set_title('A  MCIS Size vs Matching Confidence\n(NBLAST agreement as confidence proxy)',
               fontsize=10, loc='left', color='#333333')
ax_a.set_xticks(tiers)
ax_a.set_xticklabels([f'{t}%' for t in tiers])
ax_a.set_ylim(0, 120)
ax_a.set_xlim(5, 105)
ax_a.spines[['top', 'right']].set_visible(False)
ax_a.tick_params(colors='#444444')
for spine in ['left', 'bottom']:
    ax_a.spines[spine].set_color('#bbbbbb')
ax_a.yaxis.grid(True, color='#e0e0e0', linewidth=0.7, zorder=0)
ax_a.set_axisbelow(True)

# Key annotation box
key_text = (
    "Key: if the result (100% pool, N=right end) were driven\n"
    "by low-confidence matches, N would drop sharply\n"
    "when restricted to high-confidence triplets only.\n"
    "Stable or increasing N → result is confidence-robust."
)
ax_a.text(
    0.03, 0.05, key_text,
    transform=ax_a.transAxes,
    fontsize=7.5, color='#555555',
    va='bottom', ha='left',
    bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='#cccccc', alpha=0.9),
)

# ── Panel B: agreement distribution ─────────────────────────────────────────
ax_b.set_facecolor('#f8f9fa')
bars = ax_b.bar(
    agree_labels, agree_counts,
    color=agree_colors, width=0.55,
    edgecolor='white', linewidth=1.5,
    zorder=2,
)
for bar, pct, count in zip(bars, agree_pcts, agree_counts):
    ax_b.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 20,
        f'{count:,}\n({pct})',
        ha='center', va='bottom',
        fontsize=9.5, fontweight='bold', color='#222222',
    )

ax_b.set_ylabel('Number of triplets', fontsize=10)
ax_b.set_title('B  NBLAST Agreement Distribution\n(confidence proxy for 2,798 triplets)',
               fontsize=10, loc='left', color='#333333')
ax_b.set_ylim(0, max(agree_counts) * 1.25)
ax_b.spines[['top', 'right']].set_visible(False)
ax_b.tick_params(colors='#444444')
for spine in ['left', 'bottom']:
    ax_b.spines[spine].set_color('#bbbbbb')
ax_b.yaxis.grid(True, color='#e0e0e0', linewidth=0.7, zorder=0)
ax_b.set_axisbelow(True)

note = ('Note: "NBLAST agreement" = NBLAST top-1 automated match\n'
        'agrees with expert-curated match (fafb_match / manc_match)')
ax_b.text(
    0.5, -0.22, note,
    transform=ax_b.transAxes,
    fontsize=7, color='#777777',
    ha='center', va='top',
)

# ── Save ─────────────────────────────────────────────────────────────────────
plt.tight_layout(pad=1.5)
fig.savefig(OUT, dpi=180, bbox_inches='tight', facecolor='white')
print(f'Saved → {OUT}')
