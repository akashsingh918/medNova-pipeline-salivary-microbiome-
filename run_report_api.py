#!/usr/bin/env python3
"""
run_report_api.py
Convenience launcher for the FastAPI service.
"""
import uvicorn

if __name__ == "__main__":
    print("🚀 Starting Mednova Report API at http://0.0.0.0:8001")
    uvicorn.run("report_api:app", host="0.0.0.0", port=8001, reload=True, log_level="info")
