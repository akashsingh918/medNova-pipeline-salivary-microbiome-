#!/usr/bin/env python3
"""
integrated_sova_health_generator.py  (lean PDF from scorer JSON)

Reads patient_reports/<sample>_report_data.json and produces a compact PDF:
  professional_reports/<sample>_Comprehensive_Report.pdf

Content:
- Title & summary (CHI, category)
- Key metrics table
- Disease risks table (if any)
- Decision trace (why this classification)
- Top genera bar chart (from raw_abundances)

Dependencies: matplotlib only (PdfPages). No ReportLab required.
"""
import argparse, os, json
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def generate_pdf(sample_id: str,
                 reports_dir: str = "patient_reports",
                 out_dir: str = "professional_reports") -> str:
    os.makedirs(out_dir, exist_ok=True)
    in_path = os.path.join(reports_dir, f"{sample_id}_report_data.json")
    if not os.path.exists(in_path):
        raise FileNotFoundError(in_path)

    data = json.load(open(in_path))
    out_pdf = os.path.join(out_dir, f"{sample_id}_Comprehensive_Report.pdf")

    with PdfPages(out_pdf) as pdf:
        # Page 1: summary
        fig = plt.figure(figsize=(8.27, 11.69))
        plt.axis("off")
        y = 0.95
        plt.text(0.5, y, "Saliva Microbiome — Comprehensive Report",
                 ha="center", va="top", fontsize=18, fontweight="bold")
        y -= 0.04
        plt.text(0.5, y, f"Sample: {sample_id} | Date: {datetime.now().strftime('%b %d, %Y')}",
                 ha="center", va="top", fontsize=11)
        y -= 0.06

        chi = data.get("clinical_health_index", 0)
        cat = data.get("health_category", {}).get("label","Unknown")
        plt.text(0.08, y, f"Clinical Health Index (CHI): {chi:.2f}", fontsize=13, fontweight="bold"); y -= 0.03
        plt.text(0.08, y, f"Category: {cat}", fontsize=12); y -= 0.03

        # Key metrics
        km = data.get("key_metrics", {})
        plt.text(0.08, y, "Key Metrics", fontsize=12.5, fontweight="bold"); y -= 0.026
        for k in ["shannon_diversity","beneficial_pathogen_log_ratio","health_balance",
                  "scfa_producer_abundance","pathogen_load","coverage_lod"]:
            if k in km:
                plt.text(0.10, y, f"• {k.replace('_',' ').title()}: {km[k]}", fontsize=10.5); y -= 0.018

        # Risks
        risks = data.get("disease_risks", {})
        y -= 0.018
        plt.text(0.08, y, "Disease Risks", fontsize=12.5, fontweight="bold"); y -= 0.026
        if risks:
            for name, info in risks.items():
                level = info.get("level","")
                marker = info.get("marker","")
                plt.text(0.10, y, f"• {name}: {level} — {marker}", fontsize=10.5); y -= 0.018
        else:
            plt.text(0.10, y, "• None detected by current panels", fontsize=10.5); y -= 0.018

        # Decision trace
        trace = data.get("decision_trace", [])
        y -= 0.018
        plt.text(0.08, y, "Decision Trace", fontsize=12.5, fontweight="bold"); y -= 0.026
        if trace:
            for t in trace:
                plt.text(0.10, y, f"• {t}", fontsize=10.5); y -= 0.018
        else:
            plt.text(0.10, y, "• (no special rules fired)", fontsize=10.5); y -= 0.018

        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Page 2: Top genera bar chart
        abund = data.get("raw_abundances", {})
        if abund:
            items = sorted(abund.items(), key=lambda kv: kv[1], reverse=True)[:12]
            labels, values = zip(*items)
            fig = plt.figure(figsize=(8.27, 11.69))
            ax = fig.add_axes([0.12, 0.12, 0.76, 0.80])
            ax.barh(labels[::-1], [v*100 for v in values[::-1]])
            ax.set_xlabel("Relative abundance (%)")
            ax.set_title("Top Genera")
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    print(f"wrote {out_pdf}")
    return out_pdf

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--reports_dir", default="patient_reports")
    ap.add_argument("--out_dir", default="professional_reports")
    args = ap.parse_args()
    generate_pdf(args.sample_id, args.reports_dir, args.out_dir)

if __name__ == "__main__":
    main()
