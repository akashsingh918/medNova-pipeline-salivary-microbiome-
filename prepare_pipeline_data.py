#!/usr/bin/env python3
"""
prepare_pipeline_data.py  (ERR + LOD, resilient)

Build cohort-derived references from a BIOM-style TSV + ASV→Genus map:
  1) genus_percentiles.json : per-genus p50/p90/p97.5/p99 and LOD (5th pct of nonzero)
  2) reference_stats.json   : robust stats for shannon, bp_log_ratio, scfa, pathogen_load, health_balance

Notes
- Healthy cohort selection: pass --healthy_ids (txt file, one sample ID per line).
- LOD: per-genus 5th percentile of NONZERO abundances (detection-aware for 16S).
- Health/Disease “balance” uses Aitchison-style log-contrast of core communities.
- If salivadb_genus_weights.csv is missing, falls back to curated HEALTH_CORE/DISEASE_CORE sets.
"""
import argparse, json, os
import pandas as pd, numpy as np
from typing import Dict, Set

HEALTH_CORE: Set[str] = {
    "Streptococcus","Neisseria","Rothia","Haemophilus",
    "Granulicatella","Actinomyces","Veillonella","Gemella"
}
DISEASE_CORE: Set[str] = {
    "Fusobacterium","Porphyromonas","Tannerella","Treponema",
    "Prevotella","Peptostreptococcus","Parvimonas","Leptotrichia","Campylobacter"
}
SCFA_SET: Set[str] = {"Veillonella","Prevotella","Eubacterium","Propionibacterium","Megasphaera"}

def read_feature_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", skiprows=1, index_col=0)
    return df.div(df.sum(axis=0), axis=1).fillna(0)

def map_to_genus(rel: pd.DataFrame, asv_map_path: str) -> pd.DataFrame:
    asv_map = json.load(open(asv_map_path))
    rel = rel.copy()
    rel["genus"] = rel.index.map(asv_map).fillna("Unknown")
    return rel.groupby("genus").sum()

def subset_healthy(df: pd.DataFrame, healthy_txt: str | None) -> pd.DataFrame:
    if not healthy_txt:
        return df
    ids = {s.strip() for s in open(healthy_txt) if s.strip()}
    cols = [c for c in df.columns if c in ids]
    return df[cols] if cols else df

def compute_percentiles(genus_rel: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    for genus, row in genus_rel.iterrows():
        vals = row.values.astype(float)
        nz = vals[vals > 0]
        lod = float(np.nanpercentile(nz, 5)) if nz.size else 0.0
        p50, p90, p975, p99 = np.nanpercentile(vals, [50, 90, 97.5, 99]).tolist()
        out[genus] = {"p50": float(p50), "p90": float(p90),
                      "p97_5": float(p975), "p99": float(p99), "lod": lod}
    return out

def robust_stats_from_samples(genus_rel: pd.DataFrame, weights_csv: str | None) -> Dict[str, Dict[str, float]]:
    # Beneficial/pathogenic sets (weights optional)
    beneficial, pathogenic = set(HEALTH_CORE), set(DISEASE_CORE)
    if weights_csv and os.path.exists(weights_csv):
        w = pd.read_csv(weights_csv, index_col=0)
        beneficial = set(w[w["HealthWeight"] > 0].index) or beneficial
        pathogenic = set(w[w["HealthWeight"] < 0].index) or pathogenic

    eps = 1e-6
    vals = {"shannon": [], "bp_log_ratio": [], "scfa": [], "pathogen_load": [], "health_balance": []}

    for col in genus_rel.columns:
        a = genus_rel[col]
        p = a[a > 0].values
        shannon = float(-np.sum(p * np.log(p))) if p.size else 0.0
        benef = float(a[a.index.isin(beneficial)].sum())
        patho = float(a[a.index.isin(pathogenic)].sum())
        scfa = float(a[a.index.isin(SCFA_SET)].sum())
        bp_lr = float(np.log((benef + eps)/(patho + eps)))

        # Aitchison log-contrast of health vs disease cores
        hc = float(np.exp(np.mean(np.log(a.reindex(HEALTH_CORE, fill_value=0).replace(0, eps)))))
        dc = float(np.exp(np.mean(np.log(a.reindex(DISEASE_CORE, fill_value=0).replace(0, eps)))))
        balance = float(np.log((hc + eps)/(dc + eps)))

        vals["shannon"].append(shannon)
        vals["bp_log_ratio"].append(bp_lr)
        vals["scfa"].append(scfa)
        vals["pathogen_load"].append(patho)
        vals["health_balance"].append(balance)

    stats: Dict[str, Dict[str, float]] = {}
    for k, arr in vals.items():
        arr = np.asarray(arr, dtype=float)
        median = float(np.nanmedian(arr)) if len(arr) else 0.0
        mad = float(np.nanmedian(np.abs(arr - median))) if len(arr) else 1.0
        stats[k] = {
            "mean": float(np.nanmean(arr)) if len(arr) else 0.0,
            "std": float(np.nanstd(arr, ddof=1)) if len(arr) > 1 else 1.0,
            "median": median,
            "mad": mad if mad > 0 else 1.0,
            "p25": float(np.nanpercentile(arr, 25)) if len(arr) else 0.0,
            "p75": float(np.nanpercentile(arr, 75)) if len(arr) else 0.0
        }
    return stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feature_tsv", required=True)
    ap.add_argument("--asv_map", required=True)
    ap.add_argument("--healthy_ids", default=None)
    ap.add_argument("--weights_csv", default=None)
    ap.add_argument("--out_percentiles", default="genus_percentiles.json")
    ap.add_argument("--out_refstats", default="reference_stats.json")
    args = ap.parse_args()

    rel = read_feature_table(args.feature_tsv)
    genus_rel = map_to_genus(rel, args.asv_map)
    genus_h = subset_healthy(genus_rel, args.healthy_ids)

    gp = compute_percentiles(genus_h)
    json.dump(gp, open(args.out_percentiles, "w"), indent=2)

    ref = robust_stats_from_samples(genus_h, args.weights_csv)
    json.dump(ref, open(args.out_refstats, "w"), indent=2)

    print(f" wrote {args.out_percentiles} and {args.out_refstats} (healthy n={genus_h.shape[1]})")

if __name__ == "__main__":
    main()
