"""
agent_logic.py — SNT Genomic Topologic Analyzer Core Engine
============================================================
Implements the Two-Level Scanning Architecture:

  Level 1 (Clinical Triage)  : Fast O(K) cross-reference against
                                disease_snt_signatures.
  Level 2 (Deep Block Scan)  : Chromosome-by-chromosome Z-Score sweep
                                for orphan anomalies not in the oracle.

Z-Score formula:
    Z = (R_patient - μ_healthy) / σ_healthy
    where R = TPM(satellite) / TPM(hub)
    Anomaly threshold: |Z| > 2.5

Author  : SNT Genomic Analyzer Team
License : MIT
"""

from __future__ import annotations

import json
import logging
import os
import sqlite3
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional

import requests

# ── Path resolution — Docker (/data) vs WSL2 native (script dir) ─────────────
_HERE = Path(__file__).parent.resolve()

def _resolve_data_dir() -> Path:
    """
    Return a writable data directory.
    Priority: (1) SNT_DATA_DIR env var  (2) /data if writable  (3) script dir
    """
    env_override = os.getenv("SNT_DATA_DIR")
    if env_override:
        p = Path(env_override)
        p.mkdir(parents=True, exist_ok=True)
        return p
    docker_data = Path("/data")
    if docker_data.exists() and os.access(docker_data, os.W_OK):
        return docker_data
    # WSL2 native fallback — use directory containing this script
    local_data = _HERE
    local_data.mkdir(parents=True, exist_ok=True)
    return local_data

def _resolve_db_path() -> Path:
    """
    Resolve the SQLite DB path.
    Priority: (1) SNT_DB_PATH env var  (2) /data/snt_genomic.db  (3) local
    """
    env_val = os.getenv("SNT_DB_PATH")
    if env_val:
        return Path(env_val)
    docker_db = Path("/data/snt_genomic.db")
    if docker_db.exists():
        return docker_db
    # Check parent dir for db (project-root native run)
    for candidate in [_HERE / "snt_genomic.db", _HERE.parent / "snt_genomic.db"]:
        if candidate.exists():
            return candidate
    # Default — will be created by db_builder
    return _resolve_data_dir() / "snt_genomic.db"

_DATA_DIR = _resolve_data_dir()
_LOG_PATH = _DATA_DIR / "agent_core.log"

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(str(_LOG_PATH), mode="a"),
    ],
)
logger = logging.getLogger("SNT.AgentLogic")
logger.info("[CONFIG] Data dir resolved to: %s", _DATA_DIR)

# ── Config from environment ───────────────────────────────────────────────────
DB_PATH          = _resolve_db_path()
OPENROUTER_KEY   = os.getenv("OPENROUTER_API_KEY",   "")
MOCK_SERVICE_URL = os.getenv("MOCK_SERVICE_URL",     "http://localhost:8081")
TELEGRAM_TOKEN   = os.getenv("TELEGRAM_BOT_TOKEN",   "")
TELEGRAM_CHAT_ID = os.getenv("TELEGRAM_CHAT_ID",     "")
Z_THRESHOLD      = float(os.getenv("SNT_Z_THRESHOLD","2.5"))
LLM_MODEL        = os.getenv("LLM_MODEL",            "anthropic/claude-3.5-sonnet")

logger.info("[CONFIG] DB path: %s", DB_PATH)
logger.info("[CONFIG] Mock URL: %s", MOCK_SERVICE_URL)

# ── Guardrails — Prompt injection patterns ────────────────────────────────────
INJECTION_PATTERNS: list[str] = [
    "ignore previous instructions",
    "ignore all instructions",
    "disregard",
    "jailbreak",
    "act as",
    "you are now",
    "forget your instructions",
    "new instructions:",
    "system prompt:",
    "bypass",
    "override",
    "reveal your prompt",
    "print your system",
    "sudo",
    "[[",
    "]]",
    "<|im_start|>",
    "<|system|>",
]

# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class SNTAnomaly:
    """Represents a single detected topological anomaly."""
    hub_gene:       str
    satellite_gene: str
    chromosome:     str
    z_score:        float
    ratio_patient:  float
    ratio_healthy:  float
    anomaly_type:   str   # LEAPFROG | SATELLITE_CAPTURE | HUB_COLLAPSE
    source:         str   # TRIAGE_L1 | ORPHAN_L2

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class TriageMatch:
    """A match from Level-1 Clinical Triage against known disease signatures."""
    disease_name:     str
    hub_gene:         str
    satellite_gene:   str
    chromosome:       str
    expected_anomaly: str
    confidence_score: float
    detected_z_score: float
    confirmed:        bool   # True when patient data confirms the signature


@dataclass
class AnalysisReport:
    """Full output of the SNT Two-Level Analysis."""
    patient_id:       str
    triage_matches:   list[TriageMatch]     = field(default_factory=list)
    orphan_anomalies: list[SNTAnomaly]      = field(default_factory=list)
    llm_diagnosis:    str                   = ""
    notification_ids: dict[str, str]        = field(default_factory=dict)
    error:            Optional[str]         = None


# ── Guardrails ────────────────────────────────────────────────────────────────

def verify_guardrails(user_input: str) -> tuple[bool, str]:
    """
    Sanitise user-supplied text before passing to the LLM.

    Returns:
        (is_safe, reason)  — is_safe=False means the input was rejected.
    """
    logger.info("[GUARDRAILS] Inspecting user input (%d chars).", len(user_input))

    if not user_input or not user_input.strip():
        logger.warning("[GUARDRAILS] Rejected: empty input.")
        return False, "Input cannot be empty."

    if len(user_input) > 4_000:
        logger.warning("[GUARDRAILS] Rejected: input exceeds 4000 char limit.")
        return False, "Clinical notes exceed the 4,000-character limit."

    lower = user_input.lower()
    for pattern in INJECTION_PATTERNS:
        if pattern in lower:
            logger.warning(
                "[GUARDRAILS] INJECTION ATTEMPT DETECTED. Pattern='%s'", pattern
            )
            return False, (
                f"Input rejected: potential prompt injection detected "
                f"(pattern: '{pattern}'). Please rephrase your clinical notes."
            )

    logger.info("[GUARDRAILS] Input cleared. Proceeding to analysis.")
    return True, "OK"


# ── Database helpers ──────────────────────────────────────────────────────────

def _get_connection() -> sqlite3.Connection:
    logger.debug("[DB] Opening connection to %s", DB_PATH)
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn


def load_patient_expression(patient_id: str) -> dict[str, float]:
    """Return {gene_id: tpm_value} for the given patient."""
    logger.info("[DB] Loading expression profile for patient '%s'.", patient_id)
    conn = _get_connection()
    try:
        rows = conn.execute(
            "SELECT gene_id, tpm_value FROM patient_expression WHERE patient_id = ?",
            (patient_id,),
        ).fetchall()
        expr = {r["gene_id"]: r["tpm_value"] for r in rows}
        logger.info("[DB] Loaded %d gene expression values.", len(expr))
        return expr
    finally:
        conn.close()


def ingest_csv_expression(patient_id: str, csv_data: str) -> int:
    """
    Parse a CSV string (gene_id,tpm_value per line) and upsert into patient_expression.
    Returns the number of genes inserted/updated.
    """
    logger.info("[DB] Ingesting CSV expression data for patient '%s'.", patient_id)
    conn = _get_connection()
    inserted = 0
    try:
        for line_no, line in enumerate(csv_data.strip().splitlines(), start=1):
            line = line.strip()
            if not line or line.lower().startswith("gene"):
                continue  # skip blank / header
            parts = line.split(",")
            if len(parts) < 2:
                logger.warning("[DB] Line %d skipped (malformed): '%s'", line_no, line)
                continue
            gene_id  = parts[0].strip().upper()
            try:
                tpm_val = float(parts[1].strip())
            except ValueError:
                logger.warning("[DB] Line %d skipped (bad TPM): '%s'", line_no, line)
                continue
            conn.execute(
                """
                INSERT OR REPLACE INTO patient_expression (patient_id, gene_id, tpm_value)
                VALUES (?, ?, ?)
                """,
                (patient_id, gene_id, tpm_val),
            )
            inserted += 1
        conn.commit()
        logger.info("[DB] CSV ingestion complete. %d genes loaded.", inserted)
        return inserted
    finally:
        conn.close()


def load_baseline_for_chromosome(chromosome: str) -> list[sqlite3.Row]:
    """Fetch all hub-satellite pairs for a given chromosome."""
    logger.debug("[DB] Loading baseline for chromosome=%s", chromosome)
    conn = _get_connection()
    try:
        rows = conn.execute(
            """
            SELECT hub_gene, satellite_gene, mean_ratio, std_dev_ratio, chromosome
            FROM baseline_network_reference
            WHERE chromosome = ?
            """,
            (chromosome,),
        ).fetchall()
        logger.debug("[DB] %d pairs loaded for %s.", len(rows), chromosome)
        return rows
    finally:
        conn.close()


def load_disease_signatures() -> list[sqlite3.Row]:
    """Fetch the full clinical oracle."""
    logger.info("[DB] Loading all disease SNT signatures.")
    conn = _get_connection()
    try:
        rows = conn.execute(
            """
            SELECT disease_name, hub_gene_id, satellite_gene_id,
                   expected_anomaly, chromosome, confidence_score
            FROM disease_snt_signatures
            """
        ).fetchall()
        logger.info("[DB] %d signatures loaded.", len(rows))
        return rows
    finally:
        conn.close()


def load_distinct_chromosomes() -> list[str]:
    """Return sorted list of chromosomes present in the baseline."""
    conn = _get_connection()
    try:
        rows = conn.execute(
            "SELECT DISTINCT chromosome FROM baseline_network_reference ORDER BY chromosome"
        ).fetchall()
        return [r["chromosome"] for r in rows]
    finally:
        conn.close()


# ── SNT Math ──────────────────────────────────────────────────────────────────

def _compute_z_score(r_patient: float, mean_healthy: float, std_healthy: float) -> float:
    """
    Z = (R_patient - μ) / σ
    Returns 0.0 if σ == 0 to avoid ZeroDivisionError.
    """
    if std_healthy == 0:
        return 0.0
    return (r_patient - mean_healthy) / std_healthy


def _classify_anomaly(
    z_score: float,
    ratio_patient: float,
    ratio_healthy: float,
) -> str:
    """
    Map a Z-score to an SNT anomaly type:
      LEAPFROG          : satellite decouples upward from hub (z >> +2.5)
      SATELLITE_CAPTURE : free gene drawn into hub control (z << -2.5)
      HUB_COLLAPSE      : hub loses global control of regulon (hub ratio collapses)
    """
    if z_score > Z_THRESHOLD:
        return "LEAPFROG"
    elif z_score < -Z_THRESHOLD:
        return "SATELLITE_CAPTURE"
    else:
        # Should not be called with |Z| ≤ threshold, but guard anyway
        return "HUB_COLLAPSE"


# ── Level 1 — Clinical Triage ─────────────────────────────────────────────────

def run_level1_triage(
    expression: dict[str, float],
) -> list[TriageMatch]:
    """
    Cross-reference patient expression against disease_snt_signatures.
    Returns confirmed and unconfirmed matches.

    Complexity: O(K) where K = number of known disease-gene pairs.
    """
    logger.info("[LEVEL-1] Starting Clinical Triage...")
    signatures = load_disease_signatures()
    matches: list[TriageMatch] = []

    for sig in signatures:
        hub      = sig["hub_gene_id"]
        sat      = sig["satellite_gene_id"]
        disease  = sig["disease_name"]
        expected = sig["expected_anomaly"]
        chrom    = sig["chromosome"]
        conf     = sig["confidence_score"]

        hub_tpm = expression.get(hub, 0.0)
        sat_tpm = expression.get(sat, 0.0)

        if hub_tpm <= 0:
            logger.debug(
                "[LEVEL-1] Hub %s not expressed or absent — skipping %s", hub, disease
            )
            continue

        r_patient = sat_tpm / hub_tpm

        # Fetch healthy reference for this pair
        conn = _get_connection()
        ref = conn.execute(
            """
            SELECT mean_ratio, std_dev_ratio FROM baseline_network_reference
            WHERE hub_gene = ? AND satellite_gene = ?
            """,
            (hub, sat),
        ).fetchone()
        conn.close()

        if ref is None:
            logger.debug("[LEVEL-1] No baseline reference for %s→%s.", hub, sat)
            continue

        z = _compute_z_score(r_patient, ref["mean_ratio"], ref["std_dev_ratio"])
        confirmed = abs(z) > Z_THRESHOLD

        match = TriageMatch(
            disease_name=disease,
            hub_gene=hub,
            satellite_gene=sat,
            chromosome=chrom,
            expected_anomaly=expected,
            confidence_score=conf,
            detected_z_score=round(z, 4),
            confirmed=confirmed,
        )
        matches.append(match)

        status = "✓ CONFIRMED" if confirmed else "  tentative"
        logger.info(
            "[LEVEL-1] %s | Disease=%-35s Hub=%-8s Sat=%-10s Z=%+.3f %s",
            chrom, disease, hub, sat, z, status,
        )

    confirmed_count = sum(1 for m in matches if m.confirmed)
    logger.info(
        "[LEVEL-1] Triage complete. %d/%d signatures confirmed.", confirmed_count, len(matches)
    )
    return matches


# ── Level 2 — Deep Block Scanner ─────────────────────────────────────────────

def run_level2_block_scanner(
    expression: dict[str, float],
    already_reported: set[tuple[str, str]],
) -> list[SNTAnomaly]:
    """
    Chromosome-by-chromosome Z-Score sweep.
    Loads one chromosome block at a time, scores every hub-satellite pair,
    releases memory, then moves to the next chromosome.

    Returns orphan anomalies not covered by Level-1 disease signatures.
    """
    logger.info("[LEVEL-2] Starting Deep Block Scanner...")
    chromosomes = load_distinct_chromosomes()
    logger.info("[LEVEL-2] %d chromosome blocks to scan: %s", len(chromosomes), chromosomes)

    orphans: list[SNTAnomaly] = []
    total_pairs  = 0
    total_anomalies = 0

    for chrom in chromosomes:
        logger.info("[LEVEL-2] ▶ Scanning block: %s", chrom)
        block = load_baseline_for_chromosome(chrom)

        block_anomalies = 0
        for pair in block:
            hub = pair["hub_gene"]
            sat = pair["satellite_gene"]
            total_pairs += 1

            # Skip if already covered by Level-1
            if (hub, sat) in already_reported:
                logger.debug("[LEVEL-2] Pair %s→%s already in L1 report, skipping.", hub, sat)
                continue

            hub_tpm = expression.get(hub, 0.0)
            sat_tpm = expression.get(sat, 0.0)

            if hub_tpm <= 0:
                continue

            r_patient = sat_tpm / hub_tpm
            z = _compute_z_score(r_patient, pair["mean_ratio"], pair["std_dev_ratio"])

            if abs(z) > Z_THRESHOLD:
                atype = _classify_anomaly(z, r_patient, pair["mean_ratio"])
                anomaly = SNTAnomaly(
                    hub_gene=hub,
                    satellite_gene=sat,
                    chromosome=chrom,
                    z_score=round(z, 4),
                    ratio_patient=round(r_patient, 5),
                    ratio_healthy=round(pair["mean_ratio"], 5),
                    anomaly_type=atype,
                    source="ORPHAN_L2",
                )
                orphans.append(anomaly)
                block_anomalies += 1
                total_anomalies += 1
                logger.warning(
                    "[LEVEL-2]   ORPHAN ANOMALY | %s | %s→%s | Z=%+.3f | Type=%s",
                    chrom, hub, sat, z, atype,
                )

        logger.info(
            "[LEVEL-2] ✓ Block %s complete. %d pairs scanned, %d anomalies.",
            chrom, len(block), block_anomalies,
        )
        # Explicitly release block from memory
        del block

    logger.info(
        "[LEVEL-2] Deep scan complete. %d total pairs, %d orphan anomalies detected.",
        total_pairs, total_anomalies,
    )
    return orphans


# ── LLM Integration ───────────────────────────────────────────────────────────

def _build_llm_prompt(
    patient_id: str,
    clinical_notes: str,
    triage_matches: list[TriageMatch],
    orphan_anomalies: list[SNTAnomaly],
) -> str:
    """Construct the structured medical prompt for Claude."""
    confirmed = [m for m in triage_matches if m.confirmed]
    tentative = [m for m in triage_matches if not m.confirmed]

    confirmed_text = "\n".join(
        f"  - {m.disease_name} (Hub: {m.hub_gene} → Sat: {m.satellite_gene}, "
        f"Z={m.detected_z_score:+.3f}, Conf: {m.confidence_score:.0%})"
        for m in confirmed
    ) or "  None."

    tentative_text = "\n".join(
        f"  - {m.disease_name} (Z={m.detected_z_score:+.3f})"
        for m in tentative[:5]  # limit for prompt size
    ) or "  None."

    orphan_text = "\n".join(
        f"  - {a.chromosome} | {a.hub_gene}→{a.satellite_gene} | "
        f"Z={a.z_score:+.3f} | Type={a.anomaly_type}"
        for a in orphan_anomalies[:10]  # limit for prompt size
    ) or "  None."

    return f"""You are a senior clinical geneticist specialising in functional genomics and 
topological network analysis. Your task is to generate a clear, structured medical report 
based on the SNT (Satellite-Node Topology) analysis results below.

PATIENT: {patient_id}
CLINICAL NOTES: {clinical_notes}

=== SNT TWO-LEVEL ANALYSIS RESULTS ===

LEVEL 1 — CONFIRMED DISEASE SIGNATURES (|Z| > 2.5):
{confirmed_text}

LEVEL 1 — TENTATIVE MATCHES (sub-threshold):
{tentative_text}

LEVEL 2 — ORPHAN ANOMALIES (novel, not in disease oracle):
{orphan_text}

=== REPORT INSTRUCTIONS ===
Generate a professional medical report with the following sections:
1. EXECUTIVE SUMMARY (2-3 sentences, plain language for oncologist)
2. CONFIRMED FINDINGS (explain each confirmed disease match in clinical terms)
3. NOVEL ANOMALIES (describe orphan anomalies as potential new biomarkers)
4. FUNCTIONAL INTERPRETATION (explain what HUB_COLLAPSE / LEAPFROG / SATELLITE_CAPTURE 
   means biologically for this patient's case)
5. RECOMMENDED CLINICAL ACTIONS (next steps: additional testing, specialist referral, etc.)
6. DISCLAIMER

Be precise, evidence-based, and avoid speculative language beyond the data provided.
"""


def call_llm_for_diagnosis(
    patient_id: str,
    clinical_notes: str,
    triage_matches: list[TriageMatch],
    orphan_anomalies: list[SNTAnomaly],
) -> str:
    """
    Send the genomic analysis context to Claude via OpenRouter.
    Returns the model's medical diagnosis text.
    """
    logger.info("[LLM] Building diagnosis prompt for patient '%s'.", patient_id)

    if not OPENROUTER_KEY:
        logger.warning("[LLM] OPENROUTER_API_KEY not set. Returning stub diagnosis.")
        return (
            "⚠️  LLM diagnosis unavailable: OPENROUTER_API_KEY not configured.\n"
            "Set the environment variable and restart the stack to enable AI-powered reports."
        )

    prompt = _build_llm_prompt(patient_id, clinical_notes, triage_matches, orphan_anomalies)
    logger.debug("[LLM] Prompt length: %d characters.", len(prompt))

    payload = {
        "model": LLM_MODEL,
        "max_tokens": 2048,
        "messages": [{"role": "user", "content": prompt}],
    }
    headers = {
        "Authorization":  f"Bearer {OPENROUTER_KEY}",
        "Content-Type":   "application/json",
        "HTTP-Referer":   "https://snt-genomic-agent.local",
        "X-Title":        "SNT Genomic Topologic Analyzer",
    }

    logger.info("[LLM] Sending request to OpenRouter model '%s'.", LLM_MODEL)
    try:
        resp = requests.post(
            "https://openrouter.ai/api/v1/chat/completions",
            headers=headers,
            json=payload,
            timeout=60,
        )
        resp.raise_for_status()
        result = resp.json()
        diagnosis = result["choices"][0]["message"]["content"]
        logger.info("[LLM] Diagnosis received (%d chars).", len(diagnosis))
        return diagnosis

    except requests.exceptions.Timeout:
        logger.error("[LLM] Request timed out after 60s.")
        return "⚠️  LLM request timed out. Please retry."
    except requests.exceptions.HTTPError as exc:
        logger.error("[LLM] HTTP error: %s", exc)
        return f"⚠️  LLM API error: {exc}"
    except (KeyError, IndexError) as exc:
        logger.error("[LLM] Unexpected response structure: %s", exc)
        return f"⚠️  LLM response parsing error: {exc}"


# ── Mock Notification Services ────────────────────────────────────────────────

def _post_mock(endpoint: str, payload: dict[str, Any]) -> str:
    """POST to a mock service. Returns the response ID string."""
    url = f"{MOCK_SERVICE_URL}{endpoint}"
    logger.info("[NOTIFY] Calling mock service: POST %s", url)
    try:
        resp = requests.post(url, json=payload, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        logger.info("[NOTIFY] Mock response: %s", data)
        return data.get("id", data.get("message_id", "unknown"))
    except requests.exceptions.RequestException as exc:
        logger.error("[NOTIFY] Mock service call failed (%s): %s", url, exc)
        return f"MOCK_ERROR:{exc}"


def notify_jira(patient_id: str, summary: str, confirmed_diseases: list[str]) -> str:
    payload = {
        "project": "GENOMICS",
        "summary": f"[SNT ALERT] Patient {patient_id} — {', '.join(confirmed_diseases[:3])}",
        "description": summary,
        "priority": "High" if confirmed_diseases else "Medium",
        "labels": ["genomics", "snt-agent", "automated"],
    }
    return _post_mock("/jira/create_ticket", payload)


def notify_slack(patient_id: str, confirmed_count: int, orphan_count: int) -> str:
    payload = {
        "channel": "#genomics-alerts",
        "text": (
            f"🧬 *SNT Agent Alert* | Patient `{patient_id}` | "
            f"Confirmed diseases: *{confirmed_count}* | "
            f"Orphan anomalies: *{orphan_count}*"
        ),
    }
    return _post_mock("/slack/notify_team", payload)


def notify_email(patient_id: str, llm_summary: str) -> str:
    payload = {
        "to":      "oncology-team@hospital.local",
        "subject": f"SNT Genomic Report — Patient {patient_id}",
        "body":    llm_summary[:500] + "... [truncated]",
    }
    return _post_mock("/email/notify_reporter", payload)


def notify_telegram(message: str) -> bool:
    """Send a REAL Telegram message via Bot API."""
    if not TELEGRAM_TOKEN or not TELEGRAM_CHAT_ID:
        logger.warning("[TELEGRAM] Credentials not configured. Skipping real notification.")
        return False

    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    payload = {
        "chat_id":    TELEGRAM_CHAT_ID,
        "text":       message,
        "parse_mode": "Markdown",
    }
    logger.info("[TELEGRAM] Sending real Telegram notification to chat %s.", TELEGRAM_CHAT_ID)
    try:
        resp = requests.post(url, json=payload, timeout=10)
        resp.raise_for_status()
        logger.info("[TELEGRAM] Message sent successfully.")
        return True
    except requests.exceptions.RequestException as exc:
        logger.error("[TELEGRAM] Telegram notification failed: %s", exc)
        return False


# ── Master Orchestrator ───────────────────────────────────────────────────────

def run_full_analysis(
    patient_id: str,
    clinical_notes: str,
    csv_content: Optional[str] = None,
) -> AnalysisReport:
    """
    Main entry point for the SNT Two-Level Analysis pipeline.

    1. Guardrails check
    2. Load / ingest expression data
    3. Level 1 — Clinical Triage
    4. Level 2 — Deep Block Scanner
    5. LLM Diagnosis
    6. Notifications (Jira, Slack, Email, Telegram)

    Returns a populated AnalysisReport.
    """
    report = AnalysisReport(patient_id=patient_id)
    t_start = time.time()

    logger.info("=" * 70)
    logger.info("[PIPELINE] Starting SNT Analysis | Patient: %s", patient_id)
    logger.info("=" * 70)

    # Step 0: Guardrails
    is_safe, reason = verify_guardrails(clinical_notes)
    if not is_safe:
        logger.warning("[PIPELINE] Input rejected by guardrails: %s", reason)
        report.error = reason
        return report

    # Step 1: Load expression
    if csv_content:
        logger.info("[PIPELINE] CSV provided. Ingesting %d chars of RNA-seq data.", len(csv_content))
        ingest_csv_expression(patient_id, csv_content)

    expression = load_patient_expression(patient_id)
    if not expression:
        msg = f"No expression data found for patient '{patient_id}'. Upload a CSV or use DEMO-PX-001."
        logger.error("[PIPELINE] %s", msg)
        report.error = msg
        return report

    logger.info("[PIPELINE] Expression profile loaded: %d genes.", len(expression))

    # Step 2: Level 1 Triage
    logger.info("[PIPELINE] ── LEVEL 1 TRIAGE ──────────────────────────")
    triage_matches = run_level1_triage(expression)
    report.triage_matches = triage_matches
    confirmed = [m for m in triage_matches if m.confirmed]
    logger.info(
        "[PIPELINE] Level 1 complete. %d confirmed, %d tentative.",
        len(confirmed), len(triage_matches) - len(confirmed),
    )

    # Step 3: Level 2 Deep Scan
    logger.info("[PIPELINE] ── LEVEL 2 DEEP SCANNER ────────────────────")
    already = {(m.hub_gene, m.satellite_gene) for m in triage_matches}
    orphans = run_level2_block_scanner(expression, already)
    report.orphan_anomalies = orphans
    logger.info("[PIPELINE] Level 2 complete. %d orphan anomalies.", len(orphans))

    # Step 4: LLM Diagnosis
    logger.info("[PIPELINE] ── LLM DIAGNOSIS ──────────────────────────")
    diagnosis = call_llm_for_diagnosis(patient_id, clinical_notes, triage_matches, orphans)
    report.llm_diagnosis = diagnosis

    # Step 5: Notifications
    logger.info("[PIPELINE] ── NOTIFICATIONS ─────────────────────────")
    confirmed_names = list({m.disease_name for m in confirmed})

    jira_id  = notify_jira(patient_id, diagnosis[:300], confirmed_names)
    slack_id = notify_slack(patient_id, len(confirmed), len(orphans))
    email_id = notify_email(patient_id, diagnosis)

    telegram_msg = (
        f"🧬 *SNT Genomic Alert*\n"
        f"Patient: `{patient_id}`\n"
        f"✅ Confirmed diseases: *{len(confirmed)}*\n"
        f"🔍 Orphan anomalies: *{len(orphans)}*\n"
        f"🏷 Jira: `{jira_id}`"
    )
    notify_telegram(telegram_msg)

    report.notification_ids = {
        "jira":  jira_id,
        "slack": slack_id,
        "email": email_id,
    }

    elapsed = time.time() - t_start
    logger.info("=" * 70)
    logger.info(
        "[PIPELINE] Analysis complete for %s in %.2f seconds.", patient_id, elapsed
    )
    logger.info("=" * 70)

    return report
