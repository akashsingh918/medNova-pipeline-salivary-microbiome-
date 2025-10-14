#!/usr/bin/env python3
"""
report_api.py
Minimal FastAPI service to generate and download per-sample PDFs.
"""
import os, re, json, yaml
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from concurrent.futures import ThreadPoolExecutor
from subprocess import CalledProcessError, check_output

SAFE = re.compile(r"^[A-Za-z0-9_.-]+$")

app = FastAPI(title="Mednova Report API", version="1.1.0")
EXEC = ThreadPoolExecutor(max_workers=2)

REPORTS_DIR = "professional_reports"
PLOTS_DIR = "analysis_plots"
os.makedirs(REPORTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

app.mount("/reports", StaticFiles(directory=REPORTS_DIR), name="reports")
app.mount("/plots", StaticFiles(directory=PLOTS_DIR), name="plots")

@app.get("/", response_class=HTMLResponse)
def home():
    return HTMLResponse("<h2>Mednova Report API</h2><p>GET /api/generate/{sample_id}<br/>GET /api/download/{filename}</p>")

@app.get("/api/generate/{sample_id}")
def generate(sample_id: str):
    if not SAFE.match(sample_id):
        raise HTTPException(400, "Invalid sample_id")
    try:
        out = check_output([
            "python", "integrated_sova_health_generator.py",
            "--sample_id", sample_id
        ], text=True)
        # find produced pdf
        pdfs = [f for f in os.listdir(REPORTS_DIR) if f.startswith(sample_id) and f.endswith(".pdf")]
        if not pdfs:
            raise HTTPException(500, "PDF not produced")
        latest = max(pdfs, key=lambda f: os.path.getmtime(os.path.join(REPORTS_DIR, f)))
        return {"status":"ok", "pdf": f"/api/download/{latest}", "log": out}
    except CalledProcessError as e:
        raise HTTPException(500, f"Generation failed: {e.output[:4000]}")

@app.get("/api/download/{filename}")
def download(filename: str):
    if not SAFE.match(filename):
        raise HTTPException(400, "Invalid filename")
    path = os.path.join(REPORTS_DIR, filename)
    if not os.path.exists(path):
        raise HTTPException(404, "Not found")
    return FileResponse(path, filename=filename, media_type="application/pdf")

@app.get("/api/health")
def health():
    return {"status":"healthy"}
