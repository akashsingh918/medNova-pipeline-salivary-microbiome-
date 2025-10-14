prepare_pipeline_data.py

----What changed----
Builds cohort-derived genus percentiles (genus_percentiles.json) at p50/p90/p97.5 from your healthy subset.
Computes robust global stats (reference_stats.json) for Shannon, beneficial-to-pathogen log-ratio, SCFA abundance, and pathogen load using median/MAD fallbacks (safer than mean/SD).
Accepts --healthy_ids to scope the reference distributions; falls back gracefully if missing.
Normalizes feature table once (vectorized), maps ASV→genus, groups by genus.

----Why it matters-----
Replaces unverifiable “fixed weights” with Empirical Reference Ranges (ERR) from your own cohort.
Directly addresses the supervisor’s concern: e.g., healthy samples with high Fusobacterium no longer auto-flag—calls now require percentile gates relative to health.


create_health_weights.py

---What changed----
Cleans genus names (drops higher taxonomic ranks, ambiguous abbreviations like “T.”).
Enforces evidence threshold (--min_evidence, default 3) before a genus gets a HealthWeight.
Outputs salivadb_genus_weights.csv and a .meta.json with provenance (source file, timestamp, n genera).

---Why it matters---
Keeps literature-derived weights as soft priors (context), not sole arbiters of risk.
Improves reproducibility and traceability of how weights were built.


batch_health_scorer.py

----What changed----
Switches to ERR-gated multi-marker panels (OSCC, Red-complex/Periodontitis, VSC).
Panels only fire if genus abundances exceed healthy percentiles (e.g., p97.5 for Fusobacterium).
Supports “any/all/≥k” gates plus per-panel weighted sums.
Component scoring uses robust percentiles vs cohort stats; pathogen is inverted via 100 − percentile.
Composite → calibrated CHI using a piecewise mapping around cohort p5/p50/p95 to avoid overconfident extremes.
Emits clean per-sample JSON (key metrics, risks, raw genus abundances) to patient_reports/.

---Why it matters---
Dramatically reduces false positives where single taxa are naturally high in some healthy individuals.
CHI becomes interpretable and stable, grounded in the cohort distribution (not arbitrary cutoffs).


integrated_sova_health_generator.py

---What changed---
Generates a compact PDF per sample using the same scorer (no duplicated logic).
Shows summary table, top-genera bar plot, and a risks table (if any panels fired).
Requires only ReportLab + Matplotlib; plots and PDFs written to versioned folders.

---Why it matters---
Clinically friendly artifact that stays numerically faithful to the scorer.
Keeps reporting deterministic; avoids AI “narrative drift.”


report_api.py

---What changed---
Minimal FastAPI service:
GET /api/generate/{sample_id} shells out to the generator.
GET /api/download/{filename} serves PDFs.
Static mounts for /reports and /plots.
Filename whitelist (regex) to avoid shell/path shenanigans.
No business logic duplication—the API delegates to your CLI.

----Why it matters----
Simple automation surface for UI or scripts without risking logic divergence.
Safer by design; easy to debug (logs show the exact generator call).


run_report_api.py

----What changed----
Convenience launcher for the API (uvicorn, reload, log level).
Prints a clear startup line with the bound URL.

----Why it matters----
Zero-friction way to run the service locally or in a container; no guesswork
