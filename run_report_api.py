#!/usr/bin/env python3
"""
run_report_api.py  (uvicorn launcher)
"""
import uvicorn

if __name__ == "__main__":
    print("🚀 Report API at http://0.0.0.0:8001  (endpoints: /health, /api/generate/{id}, /api/download/{id})")
    uvicorn.run("report_api:app", host="0.0.0.0", port=8001, reload=True, log_level="info")
