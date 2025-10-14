#!/usr/bin/env python3
"""
prepare_pipeline_data.py
Builds cohort-derived references used by the scorer:
  1) genus_percentiles.json  (p50/p90/p97.5 per genus from healthy cohort)
  2) reference_stats.json    (global stats for metrics: shannon, bp_log_ratio, scfa, pathogen_load)

Inputs (defaults can be overridden by CLI):
  - feature-table.tsv (BIOM tsv: first row comment, first col = feature ids)
  - asv_to_genus_map.json (ASV/OTU -> genus)
  - healthy_ids.txt (optional; newline-delimited sample ids for "healthy" set)
Outputs:
  - genus_percentiles.json
  - reference_stats.json
"""
import argparse, json, pandas as pd, numpy as np

def read_feature_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", skiprows=1, index_col=0)
    # relative abundance per column
    return df.div(df.sum(axis=0), axis=1)

def map_to_genus(rel: pd.DataFrame, asv_map_path: str) -> pd.DataFrame:
    asv_map = json.load(open(asv_map_path))
    rel = rel.copy()
    rel["genus"] = rel.index.map(asv_map).fillna("Unknown")
    return rel.groupby("genus").sum()

def pick_healthy(df: pd.DataFrame, healthy_txt: str | None) -> pd.DataFrame:
    if not healthy_txt:
        return df
    healthy = {s.strip() for s in open(healthy_txt) if s.strip()}
    cols = [c for c in df.columns if c in healthy]
    return df[cols] if cols else df

def compute_percentiles(genus_rel: pd.DataFrame) -> dict:
    out = {}
    for genus, row in genus_rel.iterrows():
        vals = row.values
        if not np.isfinite(vals).any(): 
            continue
        p50, p90, p975 = np.nanpercentile(vals, [50, 90, 97.5]).tolist()
        out[genus] = {"p50": float(p50), "p90": float(p90), "p97_5": float(p975)}
    return out

def compute_global_stats(genus_rel: pd.DataFrame, beneficial: set, pathogenic: set, scfa_set: set) -> dict:
    stats = {"shannon": {}, "bp_log_ratio": {}, "scfa": {}, "pathogen_load": {}}
    eps = 1e-6
    vals = {"shannon": [], "bp_log_ratio": [], "scfa": [], "pathogen_load": []}

    for col in genus_rel.columns:
        abund = genus_rel[col]
        p = abund[abund > 0].values
        shannon = float(-np.sum(p * np.log(p))) if p.size else 0.0
        benef = float(abund[abund.index.isin(beneficial)].sum())
        patho = float(abund[abund.index.isin(pathogenic)].sum())
        scfa = float(abund[abund.index.isin(scfa_set)].sum())
        lr = float(np.log((benef + eps) / (patho + eps)))
        vals["shannon"].append(shannon)
        vals["bp_log_ratio"].append(lr)
        vals["scfa"].append(scfa)
        vals["pathogen_load"].append(patho)

    for k, arr in vals.items():
        arr = np.asarray(arr, dtype=float)
        stats[k] = {
            "mean": float(np.nanmean(arr)),
            "std":  float(np.nanstd(arr, ddof=1) if len(arr) > 1 else 1.0),
            "median": float(np.nanmedian(arr)),
            "mad": float(np.nanmedian(np.abs(arr - np.nanmedian(arr))) if len(arr) > 0 else 1.0)
        }
    return stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feature_tsv", default="feature-table.tsv")
    ap.add_argument("--asv_map", default="asv_to_genus_map.json")
    ap.add_argument("--healthy_ids", default=None)
    ap.add_argument("--weights_csv", default="salivadb_genus_weights.csv")
    ap.add_argument("--out_percentiles", default="genus_percentiles.json")
    ap.add_argument("--out_refstats", default="reference_stats.json")
    args = ap.parse_args()

    rel = read_feature_table(args.feature_tsv)
    genus_rel = map_to_genus(rel, args.asv_map)

    # define sets for global stats using weights
    w = pd.read_csv(args.weights_csv, index_col=0)
    beneficial = set(w[w["HealthWeight"] > 0].index)
    pathogenic = set(w[w["HealthWeight"] < 0].index)
    scfa_set = {"Veillonella", "Prevotella", "Eubacterium", "Propionibacterium", "Megasphaera"}

    # healthy subset for percentiles & global stats
    genus_h = pick_healthy(genus_rel, args.healthy_ids)

    gp = compute_percentiles(genus_h)
    json.dump(gp, open(args.out_percentiles, "w"), indent=2)

    ref = compute_global_stats(genus_h, beneficial, pathogenic, scfa_set)
    json.dump(ref, open(args.out_refstats, "w"), indent=2)

    print(f"✅ wrote {args.out_percentiles} and {args.out_refstats} (healthy n={genus_h.shape[1]})")

if __name__ == "__main__":
    main()
