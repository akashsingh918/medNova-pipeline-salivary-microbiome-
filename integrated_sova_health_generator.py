#!/usr/bin/env python3
"""
integrated_sova_health_generator.py
Generates a compact PDF for a single sample using batch scorer logic inline.
Keeps dependencies light (ReportLab + Matplotlib).
"""
import os, json, argparse, pandas as pd, numpy as np, matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.platypus import BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

from batch_health_scorer import load_feature_table_to_genus_abund, analyze_sample

def plot_top(abund: dict, out_path: str, topn: int = 15):
    top = dict(sorted(abund.items(), key=lambda x: x[1], reverse=True)[:topn])
    plt.figure(figsize=(8,5))
    plt.barh(list(top.keys())[::-1], list(top.values())[::-1])
    plt.xlabel("Relative abundance")
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight"); plt.close()

def build_pdf(sample_id: str, analysis: dict, out_pdf: str, top_img: str):
    doc = BaseDocTemplate(out_pdf, pagesize=A4,
                          leftMargin=18*mm,rightMargin=18*mm,topMargin=18*mm,bottomMargin=18*mm)
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="main")
    doc.addPageTemplates([PageTemplate(id="main", frames=frame)])
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="H", parent=styles["Heading1"], textColor=colors.HexColor("#2E86AB")))

    story = []
    story.append(Paragraph("MEDNOVA DIAGNOSTICS", styles["H"]))
    story.append(Paragraph("Saliva Microbiome Health Summary", styles["Heading2"]))
    story.append(Spacer(1, 6*mm))

    a = analysis
    rows = [
        ["Sample", sample_id],
        ["CHI", a["clinical_health_index"]],
        ["Category", a["health_category"]["label"]],
        ["Shannon", a["key_metrics"]["shannon_diversity"]],
        ["B/P log-ratio", a["key_metrics"]["beneficial_pathogen_log_ratio"]],
        ["SCFA abundance", a["key_metrics"]["scfa_producer_abundance"]],
        ["Pathogen load", a["key_metrics"]["pathogen_load"]],
    ]
    tbl = Table(rows, colWidths=[45*mm, 100*mm])
    tbl.setStyle(TableStyle([("GRID",(0,0),(-1,-1),0.3,colors.grey),
                             ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#e9f2fb"))]))
    story.append(tbl)
    story.append(Spacer(1, 6*mm))

    story.append(Paragraph("Top Genera", styles["Heading2"]))
    story.append(Image(top_img, width=160*mm, height=90*mm))

    story.append(Spacer(1, 6*mm))
    story.append(Paragraph("Disease Risk Panels", styles["Heading2"]))
    if a["disease_risks"]:
        rows = [["Panel","Level","Score","Threshold","Markers"]]
        for k,v in a["disease_risks"].items():
            rows.append([k, v["level"], v["panel_score"], v["threshold"], v["marker"]])
        rtbl = Table(rows, colWidths=[55*mm,25*mm,20*mm,22*mm,45*mm])
        rtbl.setStyle(TableStyle([("GRID",(0,0),(-1,-1),0.3,colors.grey),
                                  ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#fdebd0"))]))
        story.append(rtbl)
    else:
        story.append(Paragraph("No panels exceeded empirical thresholds.", styles["Normal"]))

    doc.build(story)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--feature_tsv", default="feature-table.tsv")
    ap.add_argument("--asv_map", default="asv_to_genus_map.json")
    ap.add_argument("--weights_csv", default="salivadb_genus_weights.csv")
    ap.add_argument("--ref_stats", default="reference_stats.json")
    ap.add_argument("--genus_percentiles", default="genus_percentiles.json")
    ap.add_argument("--reports_dir", default="professional_reports")
    ap.add_argument("--plots_dir", default="analysis_plots")
    args = ap.parse_args()

    os.makedirs(args.reports_dir, exist_ok=True); os.makedirs(args.plots_dir, exist_ok=True)
    asv_map = json.load(open(args.asv_map))
    weights = pd.read_csv(args.weights_csv, index_col=0)
    ref = json.load(open(args.ref_stats)); gp = json.load(open(args.genus_percentiles))

    genus_abund = load_feature_table_to_genus_abund(args.feature_tsv, asv_map)

    if args.sample_id not in genus_abund.columns:
        raise SystemExit(f"Sample '{args.sample_id}' not found in feature table.")

    analysis = analyze_sample(genus_abund[args.sample_id], ref, weights, gp)

    top_path = os.path.join(args.plots_dir, f"{args.sample_id}_top.png")
    plot_top(analysis["raw_abundances"], top_path)

    out_pdf = os.path.join(args.reports_dir, f"{args.sample_id}_Comprehensive_Report.pdf")
    build_pdf(args.sample_id, analysis, out_pdf, top_path)
    print(f"✅ PDF → {out_pdf}")

if __name__ == "__main__":
    main()
