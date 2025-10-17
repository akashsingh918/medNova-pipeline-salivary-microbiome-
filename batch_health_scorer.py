#!/usr/bin/env python3
"""
batch_health_scorer.py  (context + override + abstain)

- Detection-aware: per-genus LOD (from percentiles) to assess coverage.
- Symbiosis-aware: health vs disease core log-contrast (Aitchison).
- ERR-gated risk panels with critical-pathogen override (Fusobacterium).
- Dominance rule: when OSCC panel fires, cap composite to avoid “healthy-looking” CHI.
- Abstain on poor coverage (“Needs review”).

Inputs:
  feature-table.tsv, asv_to_genus_map.json,
  reference_stats.json, genus_percentiles.json,
  (optional) salivadb_genus_weights.csv

Outputs:
  patient_reports/<sample>_report_data.json
"""
from __future__ import annotations
import argparse, os, json, numpy as np, pandas as pd
from math import erf, sqrt
from typing import Dict, Set

HEALTH_CORE: Set[str] = {
    "Streptococcus","Neisseria","Rothia","Haemophilus",
    "Granulicatella","Actinomyces","Veillonella","Gemella"
}
DISEASE_CORE: Set[str] = {
    "Fusobacterium","Porphyromonas","Tannerella","Treponema",
    "Prevotella","Peptostreptococcus","Parvimonas","Leptotrichia","Campylobacter"
}
FUSO_SYNERGY: Set[str] = {"Peptostreptococcus","Parvimonas","Leptotrichia","Prevotella","Campylobacter"}
SCFA_SET: Set[str] = {"Veillonella","Prevotella","Eubacterium","Propionibacterium","Megasphaera"}

def robust_percentile(x: float, stats: Dict[str,float]) -> float:
    mad = stats.get("mad", 0) or 0
    if mad > 0:
        z = 0.6745 * (x - stats["median"]) / mad
    else:
        sd = stats.get("std", 1.0) or 1.0
        z = (x - stats.get("mean", 0.0)) / sd
    p = 50 * (1 + erf(z / sqrt(2)))
    return float(np.clip(p, 0, 100))

def safe_gmean(series: pd.Series, eps: float = 1e-6) -> float:
    x = np.log(series.replace(0, eps).values)
    return float(np.exp(x.mean()))

def load_feature_table_to_genus_abund(feature_tsv: str, asv_map: Dict[str,str]) -> pd.DataFrame:
    raw = pd.read_csv(feature_tsv, sep="\t", skiprows=1, index_col=0)
    rel = raw.div(raw.sum(axis=0), axis=1).fillna(0)
    rel["genus"] = rel.index.map(asv_map).fillna("Unknown")
    return rel.groupby("genus").sum()

def load_weights_df(weights_csv: str | None) -> pd.DataFrame:
    if weights_csv and os.path.exists(weights_csv):
        df = pd.read_csv(weights_csv, index_col=0)
        if "HealthWeight" in df.columns:
            return df[["HealthWeight"]]
    # fallback to curated cores
    data = {g: 0.5 for g in HEALTH_CORE}
    data.update({g: -0.5 for g in DISEASE_CORE})
    return pd.DataFrame.from_dict(data, orient="index", columns=["HealthWeight"])

def analyze_sample(abund: pd.Series,
                   ref_stats: Dict[str,Dict[str,float]],
                   weights_df: pd.DataFrame,
                   gp: Dict[str,Dict[str,float]]) -> Dict:
    abund = (abund.clip(lower=0) / max(1e-12, abund.sum()))
    eps = 1e-6
    decision = []

    # ---- detection coverage ----
    above_lod = 0
    total_considered = 0
    for g, v in abund.items():
        lod = gp.get(g, {}).get("lod", 0.0)
        if v > 0 or lod > 0:
            total_considered += 1
            if v >= max(lod, 0.0):
                above_lod += 1
    coverage = above_lod / max(1, total_considered)
    low_coverage = coverage < 0.4

    # ---- core metrics ----
    p = abund[abund > 0].values
    shannon = float(-np.sum(p * np.log(p))) if p.size else 0.0

    beneficial = set(weights_df[weights_df["HealthWeight"] > 0].index)
    pathogenic = set(weights_df[weights_df["HealthWeight"] < 0].index)

    benef = float(abund[abund.index.isin(beneficial)].sum())
    patho = float(abund[abund.index.isin(pathogenic)].sum())
    scfa_sum = float(abund[abund.index.isin(SCFA_SET)].sum())
    bp_lr = float(np.log((benef + eps)/(patho + eps)))

    hc = safe_gmean(abund.reindex(HEALTH_CORE, fill_value=0))
    dc = safe_gmean(abund.reindex(DISEASE_CORE, fill_value=0))
    health_balance_raw = float(np.log((hc + eps)/(dc + eps)))

    # percentilize metrics
    sh = robust_percentile(shannon, ref_stats["shannon"])
    ra = robust_percentile(bp_lr,   ref_stats["bp_log_ratio"])
    hb = robust_percentile(health_balance_raw, ref_stats["health_balance"])
    sc = robust_percentile(scfa_sum,ref_stats["scfa"])
    pa = 100 - robust_percentile(patho, ref_stats["pathogen_load"])

    # ---- risk panels (ERR + overrides) ----
    fuso = float(abund.get("Fusobacterium", 0.0))
    f_p975 = gp.get("Fusobacterium", {}).get("p97_5", 0.2)
    f_p99  = gp.get("Fusobacterium", {}).get("p99", 0.25)  # safe default if cohort lacks Fuso
    sh_p25 = ref_stats["shannon"]["p25"]

    synergy_hits = sum(float(abund.get(g,0.0)) >= gp.get(g,{}).get("p90",1.0) for g in FUSO_SYNERGY)
    low_diversity = shannon <= sh_p25
    health_suppressed = hb < 40.0  # percentile-like score

    oscc_extreme = (fuso >= max(f_p99, 0.20))                # Path C: extreme Fuso triggers alone
    oscc_high = (fuso >= f_p975) and (synergy_hits >= 1 or low_diversity or health_suppressed)

    risks = {}
    if oscc_extreme:
        risks["Oral Cancer (OSCC)"] = {
            "level": "High",
            "marker": "Fusobacterium (extreme)",
            "panel_score": round(fuso, 4),
            "threshold": float(f_p975)
        }
        decision.append("OSCC: extreme Fusobacterium (≥p99 or ≥0.20)")
    elif oscc_high:
        risks["Oral Cancer (OSCC)"] = {
            "level": "High" if synergy_hits >= 2 else "Moderate",
            "marker": "Fusobacterium + co-factors",
            "panel_score": round(fuso + 0.5*sum(abund.get(g,0.0) for g in FUSO_SYNERGY), 4),
            "threshold": float(f_p975)
        }
        decision.append("OSCC: high Fusobacterium (≥p97.5) + synergy/low-diversity/health-core suppression")

    # Periodontitis (Red complex): ≥2/3 genera ≥ p90
    red = ["Porphyromonas","Tannerella","Treponema"]
    red_hits = sum(float(abund.get(g,0.0)) >= gp.get(g,{}).get("p90",1.0) for g in red)
    if red_hits >= 2:
        risks["Periodontitis (Red Complex)"] = {
            "level": "Moderate" if red_hits == 2 else "High",
            "marker": "Porphyromonas/Tannerella/Treponema",
            "panel_score": round(sum(abund.get(g,0.0) for g in red), 4),
            "threshold": 0.02
        }
        decision.append(f"Periodontitis: {red_hits}/3 red-complex genera ≥ p90")

    # Halitosis / VSC
    if float(abund.get("Solobacterium",0.0)) >= gp.get("Solobacterium",{}).get("p90",1.0) or \
       (fuso >= f_p975 and low_diversity):
        risks["Halitosis (VSC)"] = {
            "level": "Moderate",
            "marker": "Solobacterium/Fusobacterium context",
            "panel_score": round(float(abund.get("Solobacterium",0.0)) + 0.5*fuso, 4),
            "threshold": 0.01
        }
        decision.append("Halitosis: Solobacterium ≥ p90 or high Fuso + low diversity")

    # ---- composite scoring with risk dominance ----
    weights = {"shannon":0.20,"ratio":0.25,"health_balance":0.20,"scfa":0.15,"pathogen":0.20}
    composite = (weights["shannon"]*sh + weights["ratio"]*ra +
                 weights["health_balance"]*hb + weights["scfa"]*sc +
                 weights["pathogen"]*pa)

    if "Oral Cancer (OSCC)" in risks:
        composite = min(composite, 50.0)  # cap to prevent “healthy” CHI
        decision.append("Dominance: OSCC risk capped composite at 50")

    # map composite→CHI ([-10,10]) with a gentle slope
    x, y = [35.0, 50.0, 70.0], [-4.0, 0.0, +4.0]
    chi = float(np.clip(np.interp(composite, x, y), -10, 10))

    abstain = low_coverage
    if abstain:
        decision.append(f"Abstain: low coverage ({coverage:.0%} genera ≥ LOD)")

    category = ("Needs review" if abstain else
                "Excellent" if chi >= 3 else "Good" if chi >= 1 else
                "Average" if chi >= -1 else "Below Average" if chi >= -3 else "Non-Ideal")

    return {
        "clinical_health_index": round(chi, 2),
        "health_category": {
            "label": category,
            "color": ("grey" if abstain else "green" if category in {"Excellent","Good"}
                      else "yellow" if category=="Average" else "red")
        },
        "key_metrics": {
            "shannon_diversity": round(shannon,3),
            "beneficial_pathogen_log_ratio": round(bp_lr,3),
            "health_balance": round(health_balance_raw,3),
            "scfa_producer_abundance": round(scfa_sum,4),
            "pathogen_load": round(patho,4),
            "coverage_lod": round(coverage,3)
        },
        "composite_score": round(composite,1),
        "disease_risks": risks,
        "raw_abundances": {k: float(v) for k,v in abund.items()},
        "decision_trace": decision
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feature_tsv", required=True)
    ap.add_argument("--asv_map", required=True)
    ap.add_argument("--ref_stats", default="reference_stats.json")
    ap.add_argument("--genus_percentiles", default="genus_percentiles.json")
    ap.add_argument("--weights_csv", default=None)  # optional
    ap.add_argument("--outdir", default="patient_reports")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    asv_map = json.load(open(args.asv_map))
    ref_stats = json.load(open(args.ref_stats))
    gp = json.load(open(args.genus_percentiles))
    weights = load_weights_df(args.weights_csv)

    genus_abund = load_feature_table_to_genus_abund(args.feature_tsv, asv_map)

    for sample_id in genus_abund.columns:
        res = analyze_sample(genus_abund[sample_id], ref_stats, weights, gp)
        res["sample_id"] = sample_id
        out = os.path.join(args.outdir, f"{sample_id}_report_data.json")
        json.dump(res, open(out, "w"), indent=2)
        print(f"✔ {sample_id} → {out}")

    print("✅ batch scoring complete.")

if __name__ == "__main__":
    main()
