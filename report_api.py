#!/usr/bin/env python3
"""
report_api.py  (thin FastAPI wrapper)

Endpoints:
- GET /api/generate/{sample_id}  → builds PDF from scorer JSON
- GET /api/download/{sample_id}  → returns latest PDF for the sample (if exists)
- GET /health                   → health check

Serves /reports as static files.
"""
import os, glob
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from datetime import datetime

from integrated_sova_health_generator import generate_pdf

REPORTS_DIR = "professional_reports"
PATIENT_JSON_DIR = "patient_reports"

app = FastAPI(title="Mednova Report API", version="1.0.0")
app.mount("/reports", StaticFiles(directory=REPORTS_DIR), name="reports")

@app.get("/health")
def health():
    return {"status":"ok","time":datetime.now().isoformat()}

@app.get("/api/generate/{sample_id}")
def api_generate(sample_id: str):
    json_path = os.path.join(PATIENT_JSON_DIR, f"{sample_id}_report_data.json")
    if not os.path.exists(json_path):
        raise HTTPException(status_code=404, detail=f"No scorer JSON for {sample_id}. Run batch_health_scorer.py first.")
    try:
        pdf_path = generate_pdf(sample_id, reports_dir=PATIENT_JSON_DIR, out_dir=REPORTS_DIR)
        return {"status":"ok","pdf":pdf_path}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/download/{sample_id}")
def api_download(sample_id: str):
    pats = sorted(glob.glob(os.path.join(REPORTS_DIR, f"{sample_id}_*.pdf")))
    if not pats:
        raise HTTPException(status_code=404, detail="No PDF found for sample")
    latest = max(pats, key=os.path.getmtime)
    return FileResponse(latest, media_type="application/pdf", filename=os.path.basename(latest))
