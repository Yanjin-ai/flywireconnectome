"""
Derived statistics — recompute every composition / enrichment / dimorphism
number quoted in science.md §5-§7 directly from the canonical circuit
(network_enriched.csv) and the FAFB annotation background.

Run:
    MCIS_DATA_DIR=/path/to/data python src/derived_stats.py

Writes results/derived_stats.json and prints a human-readable summary.
"""
import json
import os
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

from mcis_paths import data_dir, repo_root
REPO = Path(repo_root())
DATA_DIR = data_dir()


def fafb_background():
    a = pd.read_table(DATA_DIR + "fafb_annotations.tsv", low_memory=False)
    vc = a["super_class"].value_counts()
    return int(len(a)), vc.to_dict()


def annotation_quality(circuit_banc_ids):
    """Manual-annotation rate (non-null `manual_cluster`) for circuit vs
    non-circuit neurons within the matched-triplet pool, with Fisher exact p."""
    m = pd.read_feather(DATA_DIR + "banc_meta.feather")
    m["root_626"] = m["root_626"].astype(str)
    m = m.drop_duplicates("root_626")
    pool = m[m["fafb_match"].notna() & m["manc_match"].notna()]
    ic = pool["root_626"].isin(set(circuit_banc_ids))
    nn = pool["manual_cluster"].notna()
    a, b = int(nn[ic].sum()), int(ic.sum())
    c, d = int(nn[~ic].sum()), int((~ic).sum())
    _, p = fisher_exact([[a, b - a], [c, d - c]], alternative="greater")
    return {
        "metric": "manual_cluster non-null",
        "circuit": a, "circuit_n": b, "circuit_pct": round(100 * a / b, 1),
        "noncircuit": c, "noncircuit_n": d, "noncircuit_pct": round(100 * c / d, 1),
        "fisher_p": p,
    }


def main():
    enr = pd.read_csv(REPO / "network_enriched.csv")
    n = len(enr)
    print(f"Circuit N = {n}")

    # ── composition by super_class ──────────────────────────────────────
    comp = enr["super_class"].value_counts()
    comp_pct = {k: (int(v), round(100 * v / n, 1)) for k, v in comp.items()}
    print("\nComposition:")
    for k, (c, p) in comp_pct.items():
        print(f"  {k:<20s} {c:>3d}  ({p}%)")

    # ── neurotransmitter distribution ──────────────────────────────────
    nt = enr["neurotransmitter_predicted"].value_counts()
    nt_pct = {k: (int(v), round(100 * v / n, 1)) for k, v in nt.items()}
    print("\nNeurotransmitter:")
    for k, (c, p) in nt_pct.items():
        print(f"  {k:<16s} {c:>3d}  ({p}%)")

    # ── sexual dimorphism ──────────────────────────────────────────────
    dim = enr["sexually_dimorphic"].value_counts(dropna=False)
    iso = int(dim.get("isomorphic", 0))
    dimorphic = n - iso
    print(f"\nSexual conservation: {iso}/{n} isomorphic "
          f"({100*iso/n:.1f}%); {dimorphic} dimorphic")

    # ── cell-type enrichment vs FAFB background (Fisher exact) ──────────
    bg_total, bg = fafb_background()
    enrichment = {}
    for cls in ["descending", "ascending"]:
        circ = int(comp.get(cls, 0))
        bg_n = int(bg.get(cls, 0))
        # 2x2: [in-class in-circuit, in-class background] / [rest...]
        table = [[circ, bg_n - circ],
                 [n - circ, bg_total - bg_n - (n - circ)]]
        odds, p = fisher_exact(table, alternative="greater")
        circ_frac = circ / n
        bg_frac = bg_n / bg_total
        fold = circ_frac / bg_frac if bg_frac else float("nan")
        enrichment[cls] = {
            "circuit": circ, "circuit_pct": round(100 * circ_frac, 1),
            "background": bg_n, "background_pct": round(100 * bg_frac, 2),
            "fold": round(fold, 1), "p_value": p,
        }
        print(f"\n{cls.capitalize()} enrichment: "
              f"{circ}/{n} ({100*circ_frac:.1f}%) vs "
              f"{bg_n}/{bg_total} ({100*bg_frac:.2f}%) "
              f"→ {fold:.1f}×  (Fisher p={p:.2e})")

    annot = annotation_quality(enr["BANC"].astype(str))
    print(f"\nAnnotation quality (manual_cluster): "
          f"circuit {annot['circuit_pct']}% vs non-circuit {annot['noncircuit_pct']}% "
          f"(Fisher p={annot['fisher_p']:.3f})")

    out = {
        "N": n,
        "composition": comp_pct,
        "neurotransmitter": nt_pct,
        "sexual_conservation": {
            "isomorphic": iso, "dimorphic": dimorphic,
            "pct_isomorphic": round(100 * iso / n, 1),
        },
        "enrichment": enrichment,
        "annotation_quality": annot,
        "fafb_background_total": bg_total,
    }
    (REPO / "results").mkdir(exist_ok=True)
    with open(REPO / "results" / "derived_stats.json", "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nWrote {REPO/'results'/'derived_stats.json'}")


if __name__ == "__main__":
    main()
