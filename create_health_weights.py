#!/usr/bin/env python3
"""
create_health_weights.py
Derive evidence-weighted genus scores from a downloaded SalivaDB CSV.

Output: salivadb_genus_weights.csv with a 'HealthWeight' column in [-1, 1].
"""
import pandas as pd, numpy as np, os, re, argparse, json
from datetime import datetime

EXCLUDE_TAXA = {
    "Bacteria","Eubacteria","Proteobacteria","Firmicutes","Actinobacteria","Bacteroidetes",
    "Fusobacteria","Spirochaetes","Cyanobacteria","Tenericutes",
    "Bacteroidia","Flavobacteria","Betaproteobacteria","Gammaproteobacteria","Deltaproteobacteria",
    "Flavobacteriales","Spirochaetales","Bacteroidales","Clostridiales",
    "Unknown","Unclassified",""
}

def clean_genus(name: str) -> str:
    s = str(name).replace("Candidatus ", "").strip()
    s = re.sub(r"\[|\]|[a-z]__", "", s)
    tokens = s.split()
    if not tokens: return "Unknown"
    if re.match(r"^[A-Z]\.$", tokens[0]): return "Unknown"  # ambiguous abbreviations
    for t in tokens:
        if re.match(r"^[A-Z][a-z]+$", t):
            return t
    return "Unknown"

def build_weights(csv_path: str, out_csv: str, min_evidence: int = 3):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(csv_path)

    df = pd.read_csv(csv_path)
    required = {"Biomarker Name", "Regulation"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV missing columns: {required}")

    x = df[["Biomarker Name","Regulation"]].dropna()
    x = x[x["Regulation"].isin(["Upregulated","Downregulated"])].copy()
    x["Genus"] = x["Biomarker Name"].apply(clean_genus)

    counts = x.groupby(["Genus","Regulation"]).size().unstack(fill_value=0)
    if "Upregulated" not in counts: counts["Upregulated"] = 0
    if "Downregulated" not in counts: counts["Downregulated"] = 0
    counts.rename(columns={"Upregulated":"N_disease","Downregulated":"N_healthy"}, inplace=True)

    counts = counts[~counts.index.isin(EXCLUDE_TAXA)]
    counts = counts[~counts.index.str.match(r"^[A-Z]$")]

    total = counts["N_healthy"] + counts["N_disease"]
    confident = counts[total >= min_evidence].copy()
    if confident.empty:
        raise RuntimeError("No genera meet evidence threshold")

    confident["HealthWeight"] = ((confident["N_healthy"] - confident["N_disease"]) /
                                 (confident["N_healthy"] + confident["N_disease"])).fillna(0.0)

    confident[["HealthWeight"]].to_csv(out_csv)
    meta = {
        "source_csv": os.path.basename(csv_path),
        "generated_at": datetime.now().isoformat(),
        "min_evidence": min_evidence,
        "n_genera": int(len(confident)),
    }
    with open(out_csv.replace(".csv",".meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"✅ wrote {out_csv} (+ .meta.json) with {len(confident)} genera")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="SalivaDB CSV")
    ap.add_argument("--out", default="salivadb_genus_weights.csv")
    ap.add_argument("--min_evidence", type=int, default=3)
    args = ap.parse_args()
    build_weights(args.csv, args.out, args.min_evidence)
