"""
agent_logic.py — SNT Genomic Topologic Analyzer Core Engine
============================================================
Genomic implementation of the Shadow Node Theory (SNT) v2.5.0.

Base theory  : Shadow Node Theory · https://github.com/Inzainos/The-shadow-Node-Theory
SSRN preprint: https://ssrn.com/abstract=6418778
Zenodo       : https://doi.org/10.5281/zenodo.19446521
Corpus       : v5 · n=721 real cases · 89.3% sig · ρ=−0.678 · p=2.50×10⁻⁹⁷

Implements the Two-Level Scanning Architecture:

  Level 1 (Clinical Triage)  : Fast O(K) cross-reference against
                                disease_snt_signatures.
  Level 2 (Deep Block Scan)  : Chromosome-by-chromosome Z-Score sweep
                                for orphan anomalies not in the oracle.

Z-Score formula:
    Z = (R_patient - μ_healthy) / σ_healthy
    where R = TPM(satellite) / TPM(hub)
    Anomaly threshold: |Z| > 2.5

Author  : SNT Genomic Analyzer Team · Fractal Core Research
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

try:
    from scipy.optimize import curve_fit
    from scipy.stats import spearmanr
    SCIPY_OK = True
except ImportError:
    SCIPY_OK = False

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
# SECURITY: Never run DEBUG in production — exposes PHI in log files.
# Set SNT_LOG_LEVEL=INFO (or WARNING) in .env / docker-compose.yml.
_LOG_LEVEL = os.getenv("SNT_LOG_LEVEL", "INFO").upper()
_EFFECTIVE_LEVEL = getattr(logging, _LOG_LEVEL, logging.INFO)
logging.basicConfig(
    level=_EFFECTIVE_LEVEL,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(str(_LOG_PATH), mode="a"),
    ],
)
logger = logging.getLogger("SNT.AgentLogic")
logger.info("[CONFIG] Data dir resolved to: %s", _DATA_DIR)

def _log_pid(patient_id: str) -> str:
    """Return a short pseudonymous alias for a patient_id safe to write in logs.
    Uses first 8 hex chars of SHA-256. Irreversible — never reconstruct from logs.
    SECURITY NOTE: Never log patient_id in plaintext; always call _log_pid() first.
    """
    import hashlib
    return "PID-" + hashlib.sha256(patient_id.encode()).hexdigest()[:8]

# ── Config from environment ───────────────────────────────────────────────────
DB_PATH          = _resolve_db_path()
OPENROUTER_KEY   = os.getenv("OPENROUTER_API_KEY",   "")
MOCK_SERVICE_URL = os.getenv("MOCK_SERVICE_URL",     "http://localhost:8081")
TELEGRAM_TOKEN   = os.getenv("TELEGRAM_BOT_TOKEN",   "")
TELEGRAM_CHAT_ID = os.getenv("TELEGRAM_CHAT_ID",     "")
Z_THRESHOLD      = float(os.getenv("SNT_Z_THRESHOLD","2.5"))

# ── ACO-A: friction proxy genes (high expression = high friction = Regulated Decay)
# These genes represent active DNA repair, apoptosis, and cell cycle checkpoints.
# Their combined TPM level is used as a biological friction index F.
FRICTION_GENES: list[str] = [
    "TP53",    # genome guardian
    "BRCA1",   # DNA repair
    "BRCA2",   # DNA repair
    "MLH1",    # mismatch repair
    "ATM",     # DNA damage response
    "CHEK2",   # checkpoint kinase
    "RAD51",   # homologous recombination
    "FANCD2",  # Fanconi anemia pathway
    "RB1",     # cell cycle brake
    "PTEN",    # PI3K brake
]

# ACO-A collapse modes (friction × trigger × floor)
COLLAPSE_MODES = {
    "Regulated_Orbital_Decay":  "Friction alta — caída suave ley de potencias (R²≥0.85)",
    "Cracquelure_Decay":        "Fricción≈0 + gradual — fragmentación errática",
    "Floor_Arrested":           "Fricción≈0 + abrupto + piso residual — power law a piso",
    "Catastrophic_Cliff":       "Fricción≈0 + abrupto + sin piso — super-exponencial",
    "Logistic_Sweep":           "Magnitud acotada — curva S (reemplazo de variante)",
}
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
class ACOAResult:
    """Resultado del análisis ACO-A para un hub en HUB_COLLAPSE."""
    hub_gene:        str
    tau_points:      list[float]        # timepoints τ desde extinción funcional
    tpm_values:      list[float]        # valores TPM en cada τ
    delta:           float              # exponente de absorción Δ (NaN si no calculable)
    delta_r2:        float              # bondad de ajuste A(τ) = c·τ^Δ
    friction_index:  float              # F: promedio TPM de FRICTION_GENES normalizado
    has_floor:       bool               # ¿hay piso residual? (TPM final > 5% del pico)
    trigger:         str                # "abrupt" | "gradual"
    collapse_mode:   str                # uno de COLLAPSE_MODES
    absorber_gene:   str                # gen que absorbe la masa del hub (mayor ganancia TPM)
    absorber_delta_tpm: float           # ganancia TPM del absorbedor en τ_final - τ_0
    mode_description: str              # descripción del modo

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class AnalysisReport:
    """Full output of the SNT Two-Level Analysis."""
    patient_id:       str
    triage_matches:   list[TriageMatch]     = field(default_factory=list)
    orphan_anomalies: list[SNTAnomaly]      = field(default_factory=list)
    aco_results:      list[ACOAResult]      = field(default_factory=list)
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
    logger.info("[DB] Loading expression profile for patient %s.", _log_pid(patient_id))
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
    logger.info("[DB] Ingesting CSV expression data for patient %s.", _log_pid(patient_id))
    conn = _get_connection()
    inserted = 0
    try:
        for line_no, line in enumerate(csv_data.strip().splitlines(), start=1):
            line = line.strip()
            if not line or line.lower().startswith("gene"):
                continue  # skip blank / header
            parts = line.split(",")
            if len(parts) < 2:
                logger.warning("[DB] Line %d skipped (malformed): '%.40s…'", line_no, line[:40])
                continue
            gene_id  = parts[0].strip().upper()
            try:
                tpm_val = float(parts[1].strip())
            except ValueError:
                logger.warning("[DB] Line %d skipped (bad TPM): '%.40s…'", line_no, line[:40])
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
    aco_results: list["ACOAResult"] | None = None,
) -> str:
    """Construct the structured medical prompt for Claude — includes ACO-A frame.

    SECURITY NOTE: This function embeds PHI (patient_id + clinical_notes) into
    the prompt string. Never log the `prompt` variable directly. Log only
    len(prompt) or a hash. Controlled by SNT_LOG_LEVEL=INFO in production.
    """
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

    # Build ACO-A section
    aco_text = "  None (snapshot-only analysis)."
    if aco_results:
        lines = []
        for a in aco_results:
            delta_str = f"{a.delta:.3f} (R²={a.delta_r2:.2f})" if a.delta == a.delta else "N/A"
            lines.append(
                f"  - Hub={a.hub_gene} | Mode={a.collapse_mode} | "
                f"Δ={delta_str} | F={a.friction_index:.2f} | "
                f"Trigger={a.trigger} | Floor={'yes' if a.has_floor else 'no'} | "
                f"Absorber={a.absorber_gene} (+{a.absorber_delta_tpm:.1f} TPM)
"
                f"    → {a.mode_description}"
            )
        aco_text = "\n".join(lines)

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

ACO-A — ORBITAL COLLAPSE ANALYSIS (b⊥Δ framework):
{aco_text}

=== REPORT INSTRUCTIONS ===
Generate a professional medical report with the following sections:
1. EXECUTIVE SUMMARY (2-3 sentences, plain language for oncologist)
2. CONFIRMED FINDINGS (explain each confirmed disease match in clinical terms)
3. NOVEL ANOMALIES (describe orphan anomalies as potential new biomarkers)
4. FUNCTIONAL INTERPRETATION (explain what HUB_COLLAPSE / LEAPFROG / SATELLITE_CAPTURE 
   means biologically for this patient's case)
5. RECOMMENDED CLINICAL ACTIONS (next steps: additional testing, specialist referral, etc.)
6. ACO-A ORBITAL COLLAPSE INTERPRETATION (if ACO-A data present):
   - Explain the collapse mode detected (Regulated Decay vs Catastrophic Cliff)
   - Describe the absorber gene and its clinical significance
   - Interpret the Δ exponent if available (fast vs slow absorption)
   - Connect to friction: which biological mechanisms are providing friction (or not)
7. DISCLAIMER

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
    logger.info("[LLM] Building diagnosis prompt for patient %s.", _log_pid(patient_id))

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
        "X-Title":        "SNT Genomic Topologic Analyzer v2.5.0",
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



# ── ACO-A Engine ──────────────────────────────────────────────────────────────

def load_timeseries(patient_id: str, gene_id: str) -> tuple[list[float], list[float]]:
    """Fetch (timepoints, tpm_values) for a gene from timeseries table."""
    conn = _get_connection()
    try:
        rows = conn.execute(
            """
            SELECT timepoint, tpm_value
            FROM patient_expression_timeseries
            WHERE patient_id = ? AND gene_id = ?
            ORDER BY timepoint ASC
            """,
            (patient_id, gene_id),
        ).fetchall()
        if not rows:
            return [], []
        taus = [r["timepoint"] for r in rows]
        tpms = [r["tpm_value"]  for r in rows]
        return taus, tpms
    finally:
        conn.close()


def compute_friction_index(expression: dict[str, float]) -> float:
    """
    F = mean(TPM_friction_genes) / 100.0, clamped [0, 1].
    F > 0.5 → high friction (Regulated Decay)
    F ≤ 0.5 → low friction (Cliff / Cracquelure risk)
    """
    vals = [expression.get(g, 0.0) for g in FRICTION_GENES if expression.get(g, 0.0) > 0]
    if not vals:
        return 0.0
    return min(1.0, sum(vals) / len(vals) / 100.0)


def _power_law(tau: float, c: float, delta: float) -> float:
    """A(τ) = c · |τ|^Δ  (τ > 0 post-extinction)"""
    import math
    if tau <= 0:
        return c
    return c * (tau ** delta)


def fit_delta(taus: list[float], tpms: list[float]) -> tuple[float, float]:
    """
    Fit A(τ) = c · τ^Δ on post-extinction data (τ > 0).
    Returns (Δ, R²). Returns (float('nan'), 0.0) if insufficient data.
    """
    import math
    import numpy as np

    post = [(t, v) for t, v in zip(taus, tpms) if t > 0 and v > 0]
    if len(post) < 2:
        return float("nan"), 0.0

    ts = np.array([p[0] for p in post], dtype=float)
    vs = np.array([p[1] for p in post], dtype=float)

    if not SCIPY_OK:
        # Fallback: log-linear regression for power law
        log_t = np.log(ts)
        log_v = np.log(vs)
        delta_est = np.polyfit(log_t, log_v, 1)[0]
        ss_res = np.sum((log_v - np.polyval(np.polyfit(log_t, log_v, 1), log_t)) ** 2)
        ss_tot = np.sum((log_v - np.mean(log_v)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return float(delta_est), float(r2)

    try:
        popt, _ = curve_fit(
            lambda t, c, d: c * np.power(t, d),
            ts, vs,
            p0=[vs[0], -0.5],
            maxfev=5000,
        )
        c_fit, delta_fit = popt
        v_pred = c_fit * np.power(ts, delta_fit)
        ss_res = np.sum((vs - v_pred) ** 2)
        ss_tot = np.sum((vs - np.mean(vs)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return float(delta_fit), float(r2)
    except Exception:
        return float("nan"), 0.0


def classify_collapse_mode(
    friction_index: float,
    trigger: str,
    has_floor: bool,
    taus: list[float],
    tpms: list[float],
) -> str:
    """
    Classify collapse mode per ACO-A taxonomy:
      friction × trigger × floor → mode
    """
    if friction_index > 0.5:
        return "Regulated_Orbital_Decay"

    # Low friction branch
    if len(tpms) >= 2:
        # Check if magnitude is bounded (logistic — e.g. frequency-domain collapse)
        max_tpm = max(tpms)
        min_tpm = min(tpms)
        if max_tpm > 0 and (max_tpm - min_tpm) / max_tpm < 0.3:
            return "Logistic_Sweep"

    if trigger == "gradual":
        return "Cracquelure_Decay"
    else:
        # Abrupt
        if has_floor:
            return "Floor_Arrested"
        else:
            return "Catastrophic_Cliff"


def find_absorber(
    hub_gene: str,
    expression_t0: dict[str, float],
    expression_tfinal: dict[str, float],
) -> tuple[str, float]:
    """
    Find the gene with the largest TPM gain when the hub collapses.
    This is the absorber — the node capturing the hub's regulatory mass.
    Returns (absorber_gene, delta_tpm).
    """
    best_gene  = "unknown"
    best_delta = 0.0
    # Only look at genes that were satellites of this hub in the baseline
    conn = _get_connection()
    try:
        rows = conn.execute(
            "SELECT satellite_gene FROM baseline_network_reference WHERE hub_gene = ?",
            (hub_gene,),
        ).fetchall()
        satellites = {r["satellite_gene"] for r in rows}
    finally:
        conn.close()

    for gene in satellites:
        t0_val = expression_t0.get(gene, 0.0)
        tf_val = expression_tfinal.get(gene, 0.0)
        delta  = tf_val - t0_val
        if delta > best_delta:
            best_delta = delta
            best_gene  = gene

    return best_gene, best_delta


def run_aco_analysis(
    patient_id: str,
    confirmed_collapses: list[TriageMatch],
    expression: dict[str, float],
) -> list[ACOAResult]:
    """
    For each confirmed HUB_COLLAPSE, compute ACO-A frame:
      - Load timeseries for the hub gene
      - Fit Δ from post-extinction data
      - Classify collapse mode
      - Identify absorber gene

    Falls back to qualitative mode if only snapshot available.
    """
    logger.info("[ACO-A] Starting ACO-A analysis for %d collapsed hubs...",
                len(confirmed_collapses))
    results: list[ACOAResult] = []

    friction_index = compute_friction_index(expression)
    logger.info("[ACO-A] Friction index F=%.3f (%s)",
                friction_index,
                "HIGH" if friction_index > 0.5 else "LOW")

    for match in confirmed_collapses:
        if match.expected_anomaly != "HUB_COLLAPSE":
            continue

        hub = match.hub_gene
        logger.info("[ACO-A] Analyzing hub: %s (Disease: %s)", hub, match.disease_name)

        # Load timeseries
        taus, tpms = load_timeseries(patient_id, hub)
        has_timeseries = len(taus) >= 2

        # Determine trigger: abrupt if TPM drop > 70% between first two timepoints
        trigger = "gradual"
        if has_timeseries and len(tpms) >= 2:
            drop_pct = (tpms[0] - tpms[1]) / tpms[0] if tpms[0] > 0 else 0
            trigger  = "abrupt" if drop_pct > 0.5 else "gradual"
        elif not has_timeseries:
            # Infer from Z-score magnitude
            trigger = "abrupt" if abs(match.detected_z_score) > 4.0 else "gradual"

        # Floor: final TPM > 5% of initial?
        has_floor = False
        if has_timeseries and tpms:
            has_floor = tpms[-1] > 0.05 * tpms[0]
        else:
            # Snapshot proxy: current TPM > 5% of healthy baseline
            conn = _get_connection()
            ref = conn.execute(
                """SELECT mean_ratio FROM baseline_network_reference
                   WHERE hub_gene = ? LIMIT 1""",
                (hub,),
            ).fetchone()
            conn.close()
            if ref:
                healthy_approx = ref["mean_ratio"] * 100
                has_floor = expression.get(hub, 0.0) > 0.05 * healthy_approx

        # Fit Δ
        delta, delta_r2 = fit_delta(taus, tpms) if has_timeseries else (float("nan"), 0.0)

        # Classify mode
        mode = classify_collapse_mode(friction_index, trigger, has_floor, taus, tpms)

        # Find absorber — use snapshot for now (τ=0 vs τ_final proxy)
        expr_t0     = expression
        expr_tfinal = expression  # same snapshot; timeseries extends this
        if has_timeseries:
            # Build expression dict at τ_final from all timeseries genes
            conn = _get_connection()
            rows = conn.execute(
                """SELECT gene_id, tpm_value FROM patient_expression_timeseries
                   WHERE patient_id = ? AND timepoint = (
                       SELECT MAX(timepoint) FROM patient_expression_timeseries
                       WHERE patient_id = ?
                   )""",
                (patient_id, patient_id),
            ).fetchall()
            conn.close()
            expr_tfinal = {r["gene_id"]: r["tpm_value"] for r in rows} or expression

        absorber, absorber_dtpm = find_absorber(hub, expr_t0, expr_tfinal)

        result = ACOAResult(
            hub_gene        = hub,
            tau_points      = taus,
            tpm_values      = tpms,
            delta           = delta,
            delta_r2        = delta_r2,
            friction_index  = friction_index,
            has_floor       = has_floor,
            trigger         = trigger,
            collapse_mode   = mode,
            absorber_gene   = absorber,
            absorber_delta_tpm = absorber_dtpm,
            mode_description = COLLAPSE_MODES.get(mode, mode),
        )
        results.append(result)

        delta_str = f"{delta:.3f}" if delta == delta else "N/A (snapshot)"
        logger.info(
            "[ACO-A] Hub=%-8s Mode=%-28s F=%.2f Δ=%s R²=%.2f Absorber=%s (+%.1f TPM)",
            hub, mode, friction_index, delta_str, delta_r2, absorber, absorber_dtpm,
        )

    logger.info("[ACO-A] ACO-A analysis complete. %d hubs processed.", len(results))
    return results

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
    logger.info("[PIPELINE] Starting SNT Analysis | Patient: %s", _log_pid(patient_id))
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

    # Step 3b: ACO-A Analysis
    logger.info("[PIPELINE] ── ACO-A ORBITAL COLLAPSE ANALYSIS ────────")
    collapsed_hubs = [m for m in triage_matches if m.confirmed and m.expected_anomaly == "HUB_COLLAPSE"]
    aco_results = run_aco_analysis(active_patient if False else patient_id, collapsed_hubs, expression)
    report.aco_results = aco_results
    logger.info("[PIPELINE] ACO-A complete. %d hubs analyzed, modes: %s",
                len(aco_results),
                [r.collapse_mode for r in aco_results])

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
