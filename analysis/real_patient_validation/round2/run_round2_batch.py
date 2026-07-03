"""
run_round2_batch.py — Round 2: batch validation on 8 real TCGA tumor patients
=============================================================================
Second round of the real-patient validation. Runs the production SNT pipeline
(Level 1 triage, Level 2 orphan scanner, ACO-A collapse analysis) against 8
genuine TCGA-BRCA Primary Tumor cases, evaluated against the EMPIRICAL
healthy-tissue baseline (n=40 normal-adjacent samples) now shipping in
db_builder.BASELINE_NETWORK.

Data source : NIH GDC public API (open-access, de-identified), STAR-Counts,
              tpm_unstranded. Case barcodes + GDC file IDs in
              tumor_round2_file_ids.txt.
Panel data  : data/<CASE>_snt_panel_tpm.csv (58/59 SNT genes; PI3K excluded).

Reproduce:
    python run_round2_batch.py            # writes round2_summary.json
"""
import csv
import glob
import json
import os
import sqlite3
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
GENOMIC_AGENT_DIR = HERE.parents[2]
sys.path.insert(0, str(GENOMIC_AGENT_DIR / "genomic_database"))
sys.path.insert(0, str(GENOMIC_AGENT_DIR / "agent_core"))

os.environ.setdefault("SNT_LOG_LEVEL", "ERROR")
DB = Path(os.environ.get("SNT_DB_PATH", "/tmp/snt_round2.db"))
os.environ["SNT_DB_PATH"] = str(DB)
if DB.exists():
    DB.unlink()

import db_builder  # noqa: E402

conn = sqlite3.connect(str(DB))
db_builder.create_schema(conn)
db_builder.seed_baseline(conn)            # empirical baseline
db_builder.seed_disease_signatures(conn)
db_builder.seed_healing_patterns(conn)

csv_files = sorted((HERE / "data").glob("*_snt_panel_tpm.csv"))
for path in csv_files:
    case = path.name.replace("_snt_panel_tpm.csv", "")
    cur = conn.cursor()
    for row in csv.DictReader(open(path)):
        cur.execute(
            "INSERT OR REPLACE INTO patient_expression (patient_id, gene_id, tpm_value) VALUES (?,?,?)",
            (case, row["gene_id"], float(row["tpm_value"])),
        )
    conn.commit()
conn.close()

import agent_logic as A  # noqa: E402

results = []
for path in csv_files:
    case = path.name.replace("_snt_panel_tpm.csv", "")
    expr = A.load_patient_expression(case)
    triage = A.run_level1_triage(expr)
    confirmed = [m for m in triage if m.confirmed]
    already = {(m.hub_gene, m.satellite_gene) for m in triage}
    orphans = A.run_level2_block_scanner(expr, already)
    collapsed = [m for m in triage if m.confirmed and m.expected_anomaly == "HUB_COLLAPSE"]
    aco = A.run_aco_analysis(case, collapsed, expr)
    results.append({
        "case": case,
        "genes": len(expr),
        "confirmed": len(confirmed),
        "orphans": len(orphans),
        "aco_hubs": len(aco),
        "diseases": sorted({m.disease_name for m in confirmed}),
        "collapse_modes": sorted({r.collapse_mode for r in aco}),
    })
    print(f"{case}: genes={len(expr)} confirmed={len(confirmed)} "
          f"orphans={len(orphans)} aco={len(aco)}")

agg = {
    "round": 2,
    "n_patients": len(results),
    "data_source": "NIH GDC API (TCGA-BRCA Primary Tumor, open-access, de-identified)",
    "baseline": "empirical (n=40 TCGA-BRCA normal-adjacent)",
    "confirmed_mean": round(statistics.mean(r["confirmed"] for r in results), 2),
    "orphans_mean": round(statistics.mean(r["orphans"] for r in results), 2),
    "aco_mean": round(statistics.mean(r["aco_hubs"] for r in results), 2),
    "patients": results,
    "exceptions": 0,
}
(HERE / "round2_summary.json").write_text(json.dumps(agg, indent=2))
print(f"\nPatients: {agg['n_patients']} | confirmed mean: {agg['confirmed_mean']} | "
      f"orphans mean: {agg['orphans_mean']} | aco mean: {agg['aco_mean']}")
if DB.exists():
    DB.unlink()
