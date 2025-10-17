#!/usr/bin/env python3
"""
create_health_weights.py  (clean + provenance)

Derive genus-level HealthWeight from SalivaDB export:
  HealthWeight = (N_healthy - N_disease) / (N_healthy + N_disease)

- Cleans biomarker names → genus; drops higher ranks and ambiguous abbreviations.
- Applies --min_evidence threshold (default 3).
- Emits:
    salivadb_genus_weights.csv       (index=Genus, col=HealthWeight)
    salivadb_genus_weights.meta.json (provenance)
"""
import argparse, os, re, json, pandas as pd
from datetime import datetime

REQ_COLS = {"Biomarker Name", "Regulation"}
EXCLUDE_TAXA = {
    "Bacteria","Eubacteria","Proteobacteria","Firmicutes","Actinobacteria","Bacteroidetes",
    "Fusobacteria","Spirochaetes","Cyanobacteria","Tenericutes",
    "Bacteroidia","Flavobacteria","Betaproteobacteria","Gammaproteobacteria","Deltaproteobacteria",
    "Flavobacteriales","Spirochaetales","Bacteroidales","Clostridiales",
    "Unknown","Unclassified",""
}

def clean_genus(name: str) -> str:
    s = str(name).replace("Candidatus ", "").strip()
    tok = s.split()
    if not tok: return "Unknown"
    if re.match(r"^[A-Z]\.$", tok[0]):  # reject "T."
        return "Unknown"
    for t in tok:
        if re.match(r"^[A-Z][a-z]+$", t):
            return t
    return "Unknown"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="SalivaDB export CSV")
    ap.add_argument("--out_csv", default="salivadb_genus_weights.csv")
    ap.add_argument("--min_evidence", type=int, default=3)
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        raise FileNotFoundError(args.csv)

    df = pd.read_csv(args.csv)
    if not REQ_COLS.issubset(df.columns):
        raise ValueError(f"CSV missing columns: {REQ_COLS}")

    df = df[["Biomarker Name","Regulation"]].dropna()
    df = df[df["Regulation"].isin(["Upregulated","Downregulated"])]
    df["Genus"] = df["Biomarker Name"].apply(clean_genus)

    counts = df.groupby(["Genus","Regulation"]).size().unstack(fill_value=0)
    if "Upregulated" not in counts: counts["Upregulated"] = 0
    if "Downregulated" not in counts: counts["Downregulated"] = 0
    counts = counts.rename(columns={"Upregulated":"N_disease","Downregulated":"N_healthy"})

    counts = counts[~counts.index.isin(EXCLUDE_TAXA)]
    counts = counts[~counts.index.str.match(r"^[A-Z]$")]

    total = counts["N_healthy"] + counts["N_disease"]
    counts = counts[total >= args.min_evidence].copy()
    if counts.empty:
        raise RuntimeError("No genera met the evidence threshold.")

    counts["HealthWeight"] = (
        (counts["N_healthy"] - counts["N_disease"]) / (counts["N_healthy"] + counts["N_disease"])
    ).fillna(0.0)

    out = counts[["HealthWeight"]].sort_values("HealthWeight", ascending=False)
    out.to_csv(args.out_csv)

    meta = {
        "source_csv": os.path.abspath(args.csv),
        "generated_at": datetime.now().isoformat(),
        "min_evidence": args.min_evidence,
        "n_genera": int(out.shape[0]),
        "columns": list(out.columns)
    }
    with open(args.out_csv.replace(".csv",".meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"✅ wrote {args.out_csv} and {args.out_csv.replace('.csv','.meta.json')} (n={out.shape[0]})")

if __name__ == "__main__":
    main()
