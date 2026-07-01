"""
run_real_validation.py — Real-patient end-to-end validation
=============================================================
Runs the production SNT Genomic Topologic Analyzer pipeline (Level 1
clinical triage, Level 2 deep block scanner, ACO-A orbital collapse
analysis) against a genuine TCGA-BRCA patient case, instead of the
synthetic DEMO-PX-001 sample shipped with the agent.

Patient case : TCGA-BH-A18H (TCGA-BRCA cohort)
Source       : NIH GDC public API (open-access, de-identified),
               gene-expression quantification (STAR-Counts, tpm_unstranded)
File ID      : 744a6d3d-b666-49aa-8d26-47f34e3d1eb5
Gene panel   : 59/60 SNT genes recovered (PI3K excluded — pathway/family
               name, not a single gene symbol; its catalytic subunit
               PIK3CA is captured separately)

Builds a throwaway SQLite DB (path via SNT_DB_PATH, default /tmp), loads
the real TPM values, and calls the same agent_logic functions used in
production. LLM diagnosis and Jira/Slack/Email/Telegram notifications are
skipped — they require live mock services / API keys unrelated to what
this script validates (the SNT scanning/classification logic itself).

Run:
    SNT_DB_PATH=/tmp/snt_real_validation.db python run_real_validation.py
"""
import csv
import json
import os
import sqlite3
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
GENOMIC_AGENT_DIR = HERE.parents[1]
sys.path.insert(0, str(GENOMIC_AGENT_DIR / "genomic_database"))
sys.path.insert(0, str(GENOMIC_AGENT_DIR / "agent_core"))

DB_PATH = Path(os.environ.get("SNT_DB_PATH", "/tmp/snt_real_validation.db"))
REAL_CSV = HERE / "tcga_BH_A18H_snt_panel_tpm.csv"
PATIENT_ID = "TCGA-BH-A18H"

if DB_PATH.exists():
    DB_PATH.unlink()

import db_builder  # noqa: E402

conn = sqlite3.connect(str(DB_PATH))
db_builder.create_schema(conn)
db_builder.seed_baseline(conn)
db_builder.seed_disease_signatures(conn)
db_builder.seed_healing_patterns(conn)
# Intentionally NOT calling seed_demo_patient / seed_timeseries —
# this run validates against a real patient only, no synthetic data.

inserted = 0
with open(REAL_CSV, newline="") as f:
    reader = csv.DictReader(f)
    cur = conn.cursor()
    for row in reader:
        cur.execute(
            """
            INSERT OR REPLACE INTO patient_expression (patient_id, gene_id, tpm_value)
            VALUES (?, ?, ?)
            """,
            (PATIENT_ID, row["gene_id"], float(row["tpm_value"])),
        )
        inserted += 1
    conn.commit()
conn.close()
print(f"Inserted {inserted} real gene TPM rows for patient {PATIENT_ID}")

import agent_logic  # noqa: E402

expression = agent_logic.load_patient_expression(PATIENT_ID)
print(f"Loaded expression profile: {len(expression)} genes for {PATIENT_ID}")

triage_matches = agent_logic.run_level1_triage(expression)
confirmed = [m for m in triage_matches if m.confirmed]

already = {(m.hub_gene, m.satellite_gene) for m in triage_matches}
orphans = agent_logic.run_level2_block_scanner(expression, already)

collapsed_hubs = [m for m in triage_matches if m.confirmed and m.expected_anomaly == "HUB_COLLAPSE"]
aco_results = agent_logic.run_aco_analysis(PATIENT_ID, collapsed_hubs, expression)

summary = {
    "patient_id": PATIENT_ID,
    "data_source": "NIH GDC API (TCGA-BRCA, open-access, de-identified)",
    "gdc_file_id": "744a6d3d-b666-49aa-8d26-47f34e3d1eb5",
    "genes_loaded": len(expression),
    "genes_panel_total": 60,
    "level1_signatures_evaluated": len(triage_matches),
    "level1_confirmed": [
        {
            "disease": m.disease_name,
            "hub": m.hub_gene,
            "satellite": m.satellite_gene,
            "z_score": m.detected_z_score,
            "expected_anomaly": m.expected_anomaly,
        }
        for m in confirmed
    ],
    "level2_orphan_anomalies": [
        {
            "chromosome": o.chromosome,
            "hub": o.hub_gene,
            "satellite": o.satellite_gene,
            "z_score": o.z_score,
            "type": o.anomaly_type,
        }
        for o in orphans
    ],
    "aco_a_results": [
        {
            "hub": r.hub_gene,
            "collapse_mode": r.collapse_mode,
            "friction_index": r.friction_index,
            "delta": None if r.delta != r.delta else r.delta,
            "absorber_gene": r.absorber_gene,
        }
        for r in aco_results
    ],
    "pipeline_status": "completed_no_exceptions",
}

out_path = HERE / "real_validation_summary.json"
with open(out_path, "w") as f:
    json.dump(summary, f, indent=2)

print(f"\nLevel 1 confirmed: {len(confirmed)}/{len(triage_matches)}")
print(f"Level 2 orphan anomalies: {len(orphans)}")
print(f"ACO-A hubs processed: {len(aco_results)}")
print(f"Summary written to {out_path}")
