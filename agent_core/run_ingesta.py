"""
run_ingesta.py — SNT Batch Ingestion & Load Test Automation
===========================================================
Standalone script for:
  - Single patient analysis via CLI
  - Bulk stress-test (N patients in sequence)
  - Continuous watch mode (polls a directory for new CSVs)
  - JSON report export for CI pipelines

Usage examples:
  # Single analysis with demo patient
  python3 run_ingesta.py --patient DEMO-PX-001 --notes "52F, TNBC suspected"

  # Single analysis with CSV file
  python3 run_ingesta.py --patient PX-OMEGA-001 --csv paciente_omega_rnaseq.csv \
                         --notes "MYC amplification, TP53 loss"

  # Bulk load test (10 synthetic patients)
  python3 run_ingesta.py --bulk 10 --notes "Load test run"

  # Watch mode: auto-ingest any new CSV dropped into ./watch_dir/
  python3 run_ingesta.py --watch ./watch_dir --notes "Auto-ingestion"

  # Export JSON report
  python3 run_ingesta.py --patient DEMO-PX-001 --notes "Test" --json-out report.json

Author  : SNT Genomic Analyzer Team
License : MIT
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Optional

# ── Path resolution — works both in Docker (/data) and WSL2 native ─────────
_SCRIPT_DIR = Path(__file__).parent.resolve()

def _resolve_db_path() -> Path:
    """Try /data first, fall back to script directory."""
    docker_path = Path("/data/snt_genomic.db")
    local_path  = _SCRIPT_DIR / "snt_genomic.db"
    if docker_path.exists():
        return docker_path
    if local_path.exists():
        return local_path
    # Neither exists yet — default to local for first-run
    return local_path

def _resolve_log_path() -> Path:
    docker_log = Path("/data")
    if docker_log.exists() and os.access(docker_log, os.W_OK):
        return docker_log / "run_ingesta.log"
    return _SCRIPT_DIR / "run_ingesta.log"


# ── Environment bootstrap (must happen before agent_logic imports) ──────────
os.environ.setdefault("SNT_DB_PATH",        str(_resolve_db_path()))
os.environ.setdefault("MOCK_SERVICE_URL",   os.getenv("MOCK_SERVICE_URL", "http://localhost:8081"))
os.environ.setdefault("SNT_Z_THRESHOLD",    "2.5")

# ── Logging ──────────────────────────────────────────────────────────────────
LOG_PATH = _resolve_log_path()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(str(LOG_PATH), mode="a"),
    ],
)
logger = logging.getLogger("SNT.RunIngesta")

# ── Import agent pipeline ─────────────────────────────────────────────────────
try:
    # Support both: running from agent_core/ and from project root
    sys.path.insert(0, str(_SCRIPT_DIR / "agent_core"))
    sys.path.insert(0, str(_SCRIPT_DIR))
    from agent_logic import run_full_analysis, AnalysisReport
    from data_sanitizer import DataSanitizer
except ImportError as exc:
    logger.critical("Could not import agent modules: %s", exc)
    logger.critical("Run this script from the project root or agent_core/ directory.")
    sys.exit(1)


# ── Helpers ───────────────────────────────────────────────────────────────────

_SANITIZER = DataSanitizer()

def _load_csv(path: Path) -> Optional[str]:
    """Read a CSV file and return its raw string content."""
    if not path.exists():
        logger.error("[INGESTA] CSV not found: %s", path)
        return None
    raw = path.read_text(encoding="utf-8", errors="replace")
    logger.info("[INGESTA] CSV loaded: %s (%d bytes, %d lines)",
                path.name, len(raw), raw.count("\n"))
    return raw


def _sanitize_csv(raw: str, patient_id: str) -> Optional[str]:
    """Run ETL sanitizer. Returns clean CSV string or None on failure."""
    logger.info("[INGESTA] Running DataSanitizer for patient '%s'...", patient_id)
    clean_df, report = _SANITIZER.sanitize(raw, patient_id)
    logger.info("[INGESTA] Sanitizer: %s", report.summary())
    if clean_df.empty:
        logger.error("[INGESTA] Sanitizer returned empty dataset. Aborting.")
        return None
    return DataSanitizer.to_csv_string(clean_df)


def _print_report(report: AnalysisReport) -> None:
    """Pretty-print the analysis report to stdout."""
    confirmed = [m for m in report.triage_matches if m.confirmed]
    print("\n" + "=" * 70)
    print(f"  SNT ANALYSIS REPORT — Patient: {report.patient_id}")
    print("=" * 70)

    if report.error:
        print(f"  ❌ ERROR: {report.error}")
        return

    print(f"\n  ⚡ LEVEL 1 — CONFIRMED DISEASES ({len(confirmed)})")
    for m in sorted(confirmed, key=lambda x: abs(x.detected_z_score), reverse=True):
        print(f"     • {m.disease_name:<35} Hub={m.hub_gene:<8} Sat={m.satellite_gene:<10} "
              f"Z={m.detected_z_score:+.3f}  [{m.expected_anomaly}]")

    tentative = [m for m in report.triage_matches if not m.confirmed]
    if tentative:
        print(f"\n  ⚪ LEVEL 1 — TENTATIVE ({len(tentative)})")
        for m in tentative[:5]:
            print(f"     • {m.disease_name:<35} Z={m.detected_z_score:+.3f}")

    print(f"\n  🔬 LEVEL 2 — ORPHAN ANOMALIES ({len(report.orphan_anomalies)})")
    for a in sorted(report.orphan_anomalies, key=lambda x: abs(x.z_score), reverse=True):
        print(f"     • {a.chromosome:<6} {a.hub_gene:<8}→{a.satellite_gene:<10} "
              f"Z={a.z_score:+.3f}  [{a.anomaly_type}]  ⚠️  NOVEL")

    if not report.orphan_anomalies:
        print("     (none — profile matches known disease patterns only)")

    print(f"\n  📨 NOTIFICATIONS")
    for svc, nid in report.notification_ids.items():
        print(f"     • {svc.upper():<6} → {nid}")

    print(f"\n  🤖 AI DIAGNOSIS (first 400 chars)")
    snippet = report.llm_diagnosis[:400].replace("\n", " ")
    print(f"     {snippet}{'...' if len(report.llm_diagnosis) > 400 else ''}")
    print("=" * 70 + "\n")


def _export_json(report: AnalysisReport, out_path: Path) -> None:
    """Serialise the report to JSON for CI / downstream consumers."""
    data = {
        "patient_id":        report.patient_id,
        "error":             report.error,
        "confirmed_diseases": [
            {
                "disease":    m.disease_name,
                "hub":        m.hub_gene,
                "satellite":  m.satellite_gene,
                "z_score":    m.detected_z_score,
                "anomaly":    m.expected_anomaly,
                "confidence": m.confidence_score,
            }
            for m in report.triage_matches if m.confirmed
        ],
        "orphan_anomalies": [
            {
                "chromosome": a.chromosome,
                "hub":        a.hub_gene,
                "satellite":  a.satellite_gene,
                "z_score":    a.z_score,
                "type":       a.anomaly_type,
            }
            for a in report.orphan_anomalies
        ],
        "llm_diagnosis":    report.llm_diagnosis,
        "notification_ids": report.notification_ids,
    }
    out_path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("[INGESTA] JSON report saved → %s", out_path)
    print(f"  ✅ JSON report exported: {out_path}")


# ── Analysis runners ──────────────────────────────────────────────────────────

def run_single(
    patient_id: str,
    clinical_notes: str,
    csv_path: Optional[Path] = None,
    json_out: Optional[Path] = None,
) -> AnalysisReport:
    """Run one complete SNT analysis and print the result."""
    logger.info("[INGESTA] ═══ START: single analysis | patient=%s ═══", patient_id)
    t0 = time.time()

    csv_content: Optional[str] = None
    if csv_path:
        raw = _load_csv(csv_path)
        if raw:
            csv_content = _sanitize_csv(raw, patient_id)

    report = run_full_analysis(
        patient_id=patient_id,
        clinical_notes=clinical_notes,
        csv_content=csv_content,
    )

    elapsed = time.time() - t0
    logger.info("[INGESTA] ═══ DONE: %s | %.2fs ═══", patient_id, elapsed)

    _print_report(report)
    print(f"  ⏱  Analysis completed in {elapsed:.2f}s")

    if json_out:
        _export_json(report, json_out)

    return report


def run_bulk(
    n: int,
    clinical_notes: str,
    json_out: Optional[Path] = None,
) -> list[AnalysisReport]:
    """
    Run N sequential analyses using the demo patient to stress-test
    the pipeline under repeated load.
    """
    logger.info("[INGESTA] ═══ BULK LOAD TEST: %d patients ═══", n)
    reports: list[AnalysisReport] = []
    timings: list[float] = []

    for i in range(1, n + 1):
        patient_id = f"LOAD-TEST-{i:04d}"
        print(f"\n[{i}/{n}] Running {patient_id}...")
        t0 = time.time()
        report = run_full_analysis(
            patient_id="DEMO-PX-001",   # reuse demo data for load test
            clinical_notes=f"[LOAD TEST {i}/{n}] {clinical_notes}",
            csv_content=None,
        )
        elapsed = time.time() - t0
        timings.append(elapsed)
        reports.append(report)

        confirmed_count = sum(1 for m in report.triage_matches if m.confirmed)
        status = "✅" if not report.error else "❌"
        print(f"  {status} {patient_id} | Confirmed={confirmed_count} | "
              f"Orphans={len(report.orphan_anomalies)} | {elapsed:.2f}s")

    avg  = sum(timings) / len(timings)
    total = sum(timings)
    print(f"\n{'='*50}")
    print(f"  BULK TEST COMPLETE")
    print(f"  Patients : {n}")
    print(f"  Total    : {total:.2f}s")
    print(f"  Average  : {avg:.2f}s/patient")
    print(f"  Min/Max  : {min(timings):.2f}s / {max(timings):.2f}s")
    print(f"{'='*50}\n")

    logger.info(
        "[INGESTA] Bulk complete. n=%d avg=%.2fs total=%.2fs", n, avg, total
    )

    if json_out:
        bulk_data = {
            "mode":          "bulk_load_test",
            "n_patients":    n,
            "avg_seconds":   round(avg, 3),
            "total_seconds": round(total, 3),
            "results": [
                {
                    "patient_id":        r.patient_id,
                    "confirmed_diseases": sum(1 for m in r.triage_matches if m.confirmed),
                    "orphan_anomalies":   len(r.orphan_anomalies),
                    "error":              r.error,
                }
                for r in reports
            ],
        }
        json_out.write_text(json.dumps(bulk_data, indent=2), encoding="utf-8")
        print(f"  ✅ Bulk report exported: {json_out}")

    return reports


def run_watch(watch_dir: Path, clinical_notes: str) -> None:
    """
    Watch a directory for new CSV files and auto-ingest them.
    Rename processed files to .done to avoid double-processing.
    Press Ctrl+C to stop.
    """
    watch_dir.mkdir(parents=True, exist_ok=True)
    logger.info("[INGESTA] Watch mode active. Directory: %s", watch_dir)
    print(f"\n👀 Watching {watch_dir} for new CSV files... (Ctrl+C to stop)\n")

    processed: set[str] = set()

    try:
        while True:
            csv_files = sorted(watch_dir.glob("*.csv"))
            for csv_file in csv_files:
                if csv_file.name in processed:
                    continue

                patient_id = csv_file.stem.upper().replace(" ", "-")
                logger.info("[INGESTA] Watch: new file detected → %s", csv_file.name)
                print(f"  📥 New file: {csv_file.name} → patient_id={patient_id}")

                report = run_single(
                    patient_id=patient_id,
                    clinical_notes=clinical_notes,
                    csv_path=csv_file,
                )

                # Mark as processed — rename to .done
                done_path = csv_file.with_suffix(".done")
                csv_file.rename(done_path)
                processed.add(csv_file.name)
                logger.info("[INGESTA] Watch: processed → %s", done_path.name)

            time.sleep(3)  # poll every 3 seconds

    except KeyboardInterrupt:
        print("\n\n⏹  Watch mode stopped.")
        logger.info("[INGESTA] Watch mode stopped by user.")


# ── CLI ───────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="SNT Genomic Analyzer — Batch Ingestion CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--bulk",  type=int, metavar="N",
                      help="Bulk load test: run N sequential analyses")
    mode.add_argument("--watch", type=Path, metavar="DIR",
                      help="Watch directory for new CSVs (auto-ingest)")

    parser.add_argument("--patient", type=str, default="DEMO-PX-001",
                        help="Patient ID (default: DEMO-PX-001)")
    parser.add_argument("--csv",     type=Path, metavar="FILE",
                        help="Path to RNA-seq CSV file")
    parser.add_argument("--notes",   type=str,
                        default="Automated SNT analysis run via run_ingesta.py",
                        help="Clinical notes string")
    parser.add_argument("--json-out", type=Path, metavar="FILE",
                        help="Export analysis report as JSON")
    parser.add_argument("--z-threshold", type=float, default=2.5,
                        help="Z-Score anomaly threshold (default: 2.5)")
    parser.add_argument(
        "--mode", choices=["ratio", "gene"], default="ratio",
        help=(
            "Z-score mode: 'ratio' = clinical individual "
            "(TPM_sat/TPM_hub vs healthy baseline, Level-1/Level-2); "
            "'gene' = population TCGA (TPM_gene vs cohort baseline, "
            "5-Event Wall). Default: ratio."
        ),
    )

    return parser


def main() -> None:
    parser = _build_parser()
    args   = parser.parse_args()

    # Apply threshold override
    os.environ["SNT_Z_THRESHOLD"] = str(args.z_threshold)

    logger.info("=" * 60)
    logger.info("SNT run_ingesta.py starting")
    logger.info("DB path     : %s", os.environ.get("SNT_DB_PATH"))
    logger.info("Mock URL    : %s", os.environ.get("MOCK_SERVICE_URL"))
    logger.info("Z threshold : %s", os.environ.get("SNT_Z_THRESHOLD"))
    logger.info("=" * 60)

    if args.bulk:
        run_bulk(
            n=args.bulk,
            clinical_notes=args.notes,
            json_out=args.json_out,
        )

    elif args.watch:
        run_watch(
            watch_dir=args.watch,
            clinical_notes=args.notes,
        )

    else:
        run_single(
            patient_id=args.patient,
            clinical_notes=args.notes,
            csv_path=args.csv,
            json_out=args.json_out,
        )


if __name__ == "__main__":
    main()
