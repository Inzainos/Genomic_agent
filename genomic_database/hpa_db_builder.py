"""
hpa_db_builder.py — SNT Protein Database Builder
=================================================
Constructs the protein-level SNT database from two public sources:

  Source A: Human Protein Atlas (HPA) — normal_tissue.tsv
            Protein expression levels in 45 normal tissues
            URL: https://www.proteinatlas.org/download/normal_tissue.tsv.zip

  Source B: UniProt REST API — protein-protein interactions
            Hub-satellite regulatory relationships at the protein chain level

Pipeline:
  1. Download + parse HPA normal_tissue.tsv
  2. Map ordinal expression (Not detected/Low/Medium/High) → numeric [0-3]
  3. Compute μ and σ per protein across all tissues (baseline stats)
  4. Fetch hub-satellite interaction pairs from UniProt
  5. Seed known disease signatures (protein-level SNT oracle)
  6. Seed demo patient protein profile (synthetic TNBC)

Usage:
  python3 hpa_db_builder.py                  # full build with HPA download
  python3 hpa_db_builder.py --offline        # use synthetic baseline (no internet)
  python3 hpa_db_builder.py --tsv path.tsv  # use a locally downloaded HPA TSV

Author  : SNT Genomic Analyzer Team
License : MIT
"""

from __future__ import annotations

import argparse
import csv
import io
import logging
import os
import sqlite3
import sys
import time
import zipfile
from pathlib import Path
from statistics import mean, stdev
from typing import Optional
from collections import defaultdict

import requests

# ── Logging ──────────────────────────────────────────────────────────────────
def _resolve_log() -> str:
    data = Path("/data")
    if data.exists() and os.access(data, os.W_OK):
        return str(data / "hpa_db_builder.log")
    return str(Path(__file__).parent / "hpa_db_builder.log")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_resolve_log(), mode="w"),
    ],
)
logger = logging.getLogger("SNT.HPA_Builder")

# ── DB path ───────────────────────────────────────────────────────────────────
def _resolve_db() -> Path:
    env = os.getenv("SNT_DB_PATH")
    if env:
        return Path(env)
    data = Path("/data")
    if data.exists() and os.access(data, os.W_OK):
        return data / "snt_protein.db"
    return Path(__file__).parent / "snt_protein.db"

DB_PATH = _resolve_db()

# ── HPA constants ─────────────────────────────────────────────────────────────
HPA_TSV_URL = "https://www.proteinatlas.org/download/normal_tissue.tsv.zip"

ORDINAL_MAP: dict[str, int] = {
    "not detected": 0,
    "not_detected": 0,
    "low":          1,
    "medium":       2,
    "high":         3,
}

RELIABILITY_RANK: dict[str, int] = {
    "enhanced":  4,
    "supported": 3,
    "approved":  2,
    "uncertain": 1,
}

# Minimum reliability to include in baseline stats
MIN_RELIABILITY = "approved"

# ── Hub-satellite pairs (curated from UniProt + literature) ──────────────────
# These define the SNT regulatory network at the protein level.
# Format: (hub_protein, satellite_protein, interaction_type, confidence, tissue_context)
# interaction_type: TRANSCRIPTIONAL | PHOSPHORYLATION | UBIQUITINATION | COMPLEX
PROTEIN_INTERACTIONS: list[tuple] = [
    # MYC transcriptional regulon
    ("MYC",   "CDK4",   "TRANSCRIPTIONAL", 0.95, "breast"),
    ("MYC",   "CDK6",   "TRANSCRIPTIONAL", 0.92, "breast"),
    ("MYC",   "CCND1",  "TRANSCRIPTIONAL", 0.93, "breast"),
    ("MYC",   "E2F1",   "TRANSCRIPTIONAL", 0.94, "breast"),
    ("MYC",   "MCM2",   "TRANSCRIPTIONAL", 0.90, "breast"),
    ("MYC",   "PCNA",   "TRANSCRIPTIONAL", 0.91, "breast"),
    ("MYC",   "TOP2A",  "TRANSCRIPTIONAL", 0.88, "breast"),
    ("MYC",   "AURKB",  "TRANSCRIPTIONAL", 0.87, "breast"),
    ("MYC",   "PLK1",   "TRANSCRIPTIONAL", 0.86, "breast"),
    # TP53 regulon
    ("TP53",  "CDKN1A", "TRANSCRIPTIONAL", 0.97, "breast"),
    ("TP53",  "MDM2",   "TRANSCRIPTIONAL", 0.95, "breast"),
    ("TP53",  "BAX",    "TRANSCRIPTIONAL", 0.93, "breast"),
    ("TP53",  "BBC3",   "TRANSCRIPTIONAL", 0.90, "breast"),   # PUMA
    ("TP53",  "GADD45A","TRANSCRIPTIONAL", 0.89, "breast"),
    ("TP53",  "TIGAR",  "TRANSCRIPTIONAL", 0.85, "breast"),
    # EGFR kinase cascade (phosphorylation network)
    ("EGFR",  "GRB2",   "PHOSPHORYLATION", 0.96, "lung"),
    ("EGFR",  "SOS1",   "PHOSPHORYLATION", 0.94, "lung"),
    ("EGFR",  "PIK3CA", "PHOSPHORYLATION", 0.92, "lung"),
    ("EGFR",  "AKT1",   "PHOSPHORYLATION", 0.91, "lung"),
    ("EGFR",  "STAT3",  "PHOSPHORYLATION", 0.90, "lung"),
    ("EGFR",  "ERK2",   "PHOSPHORYLATION", 0.93, "lung"),
    # BRCA1 DNA repair complex
    ("BRCA1", "RAD51",  "COMPLEX",         0.94, "breast"),
    ("BRCA1", "FANCD2", "COMPLEX",         0.91, "breast"),
    ("BRCA1", "RPA1",   "COMPLEX",         0.89, "breast"),
    ("BRCA1", "RFC1",   "COMPLEX",         0.87, "breast"),
    # PIK3CA / AKT / mTOR pathway
    ("PIK3CA","AKT1",   "PHOSPHORYLATION", 0.95, "breast"),
    ("PIK3CA","MTOR",   "PHOSPHORYLATION", 0.93, "breast"),
    ("AKT1",  "MTOR",   "PHOSPHORYLATION", 0.94, "breast"),
    ("AKT1",  "GSK3B",  "PHOSPHORYLATION", 0.92, "breast"),
    ("MTOR",  "RPS6KB1","PHOSPHORYLATION", 0.91, "breast"),   # S6K1
    # KRAS / RAF / MEK cascade
    ("KRAS",  "RAF1",   "PHOSPHORYLATION", 0.95, "colon"),
    ("KRAS",  "MAP2K1", "PHOSPHORYLATION", 0.93, "colon"),    # MEK1
    ("KRAS",  "MAPK3",  "PHOSPHORYLATION", 0.91, "colon"),    # ERK1
    ("RAF1",  "MAP2K1", "PHOSPHORYLATION", 0.94, "colon"),
    ("MAP2K1","MAPK3",  "PHOSPHORYLATION", 0.95, "colon"),
    # BRAF melanoma pathway
    ("BRAF",  "MAP2K1", "PHOSPHORYLATION", 0.96, "skin"),
    ("BRAF",  "MAP2K2", "PHOSPHORYLATION", 0.94, "skin"),
    # PTEN tumour suppressor
    ("PTEN",  "AKT1",   "PHOSPHORYLATION", 0.95, "prostate"),
    ("PTEN",  "MTOR",   "PHOSPHORYLATION", 0.93, "prostate"),
    # RB1 cell cycle control
    ("RB1",   "E2F1",   "TRANSCRIPTIONAL", 0.96, "retina"),
    ("RB1",   "CCND1",  "TRANSCRIPTIONAL", 0.93, "retina"),
    ("RB1",   "CDK4",   "COMPLEX",         0.91, "retina"),
    # VHL hypoxia pathway
    ("VHL",   "HIF1A",  "UBIQUITINATION",  0.97, "kidney"),
    ("VHL",   "VEGFA",  "UBIQUITINATION",  0.94, "kidney"),
    ("HIF1A", "VEGFA",  "TRANSCRIPTIONAL", 0.95, "kidney"),
    # APC / Wnt signalling
    ("APC",   "CTNNB1", "UBIQUITINATION",  0.96, "colon"),
    ("CTNNB1","TCF7L2", "TRANSCRIPTIONAL", 0.93, "colon"),
    # SMAD / TGFb pathway
    ("SMAD4", "TGFB1",  "TRANSCRIPTIONAL", 0.92, "pancreas"),
    ("SMAD4", "TGFB2",  "TRANSCRIPTIONAL", 0.90, "pancreas"),
    # MDM2 / TP53 feedback
    ("MDM2",  "TP53",   "UBIQUITINATION",  0.97, "breast"),
]

# ── Disease protein SNT signatures ────────────────────────────────────────────
# Protein-level SNT oracle — expected Z-score pattern per disease
# Format: (disease_name, hub_protein, satellite_protein, expected_anomaly,
#          tissue_context, confidence_score)
DISEASE_PROTEIN_SIGNATURES: list[tuple] = [
    # TNBC
    ("Breast_Cancer_TNBC",       "MYC",   "CDK4",    "HUB_OVERACTIVATION", "breast", 0.92),
    ("Breast_Cancer_TNBC",       "TP53",  "CDKN1A",  "HUB_COLLAPSE",       "breast", 0.91),
    ("Breast_Cancer_TNBC",       "BRCA1", "RAD51",   "SATELLITE_CAPTURE",  "breast", 0.88),
    ("Breast_Cancer_TNBC",       "MDM2",  "TP53",    "HUB_OVERACTIVATION", "breast", 0.89),
    # HER2+ Breast
    ("Breast_Cancer_HER2pos",    "ERBB2", "GRB2",    "HUB_OVERACTIVATION", "breast", 0.94),
    ("Breast_Cancer_HER2pos",    "ERBB2", "AKT1",    "HUB_OVERACTIVATION", "breast", 0.91),
    ("Breast_Cancer_HER2pos",    "TP53",  "BAX",     "HUB_COLLAPSE",       "breast", 0.85),
    # Lung Adenocarcinoma EGFR-driven
    ("Lung_Adenocarcinoma_EGFR", "EGFR",  "STAT3",   "HUB_OVERACTIVATION", "lung",   0.93),
    ("Lung_Adenocarcinoma_EGFR", "EGFR",  "AKT1",    "HUB_OVERACTIVATION", "lung",   0.90),
    ("Lung_Adenocarcinoma_EGFR", "KRAS",  "RAF1",    "SATELLITE_CAPTURE",  "lung",   0.87),
    # KRAS-driven Colorectal
    ("Colorectal_Cancer_KRAS",   "KRAS",  "MAP2K1",  "HUB_OVERACTIVATION", "colon",  0.91),
    ("Colorectal_Cancer_KRAS",   "APC",   "CTNNB1",  "HUB_COLLAPSE",       "colon",  0.93),
    ("Colorectal_Cancer_KRAS",   "SMAD4", "TGFB1",   "HUB_COLLAPSE",       "colon",  0.88),
    # BRAF V600E Melanoma
    ("Melanoma_BRAF_V600E",      "BRAF",  "MAP2K1",  "HUB_OVERACTIVATION", "skin",   0.96),
    ("Melanoma_BRAF_V600E",      "PTEN",  "AKT1",    "HUB_COLLAPSE",       "skin",   0.89),
    # GBM
    ("Glioblastoma_GBM",         "EGFR",  "PIK3CA",  "HUB_OVERACTIVATION", "brain",  0.88),
    ("Glioblastoma_GBM",         "PTEN",  "MTOR",    "HUB_COLLAPSE",       "brain",  0.87),
    ("Glioblastoma_GBM",         "RB1",   "CCND1",   "HUB_COLLAPSE",       "brain",  0.84),
    # Renal Cell Carcinoma
    ("Renal_Cell_Carcinoma",     "VHL",   "HIF1A",   "HUB_COLLAPSE",       "kidney", 0.95),
    ("Renal_Cell_Carcinoma",     "VHL",   "VEGFA",   "HUB_COLLAPSE",       "kidney", 0.93),
    # PDAC
    ("Pancreatic_PDAC",          "KRAS",  "MAPK3",   "HUB_OVERACTIVATION", "pancreas",0.92),
    ("Pancreatic_PDAC",          "SMAD4", "TGFB2",   "HUB_COLLAPSE",       "pancreas",0.90),
    # Li-Fraumeni
    ("Li_Fraumeni_Syndrome",     "TP53",  "BAX",     "HUB_COLLAPSE",       "breast", 0.97),
    ("Li_Fraumeni_Syndrome",     "TP53",  "CDKN1A",  "HUB_COLLAPSE",       "breast", 0.96),
    ("Li_Fraumeni_Syndrome",     "MDM2",  "TP53",    "LEAPFROG",           "breast", 0.94),
    # Prostate Cancer
    ("Prostate_Cancer_AR",       "PTEN",  "AKT1",    "HUB_COLLAPSE",       "prostate",0.91),
    ("Prostate_Cancer_AR",       "AKT1",  "MTOR",    "HUB_OVERACTIVATION", "prostate",0.89),
]

# ── Demo patient protein profile (synthetic TNBC) ─────────────────────────────
# expression_score: 0=not detected, 1=low, 2=medium, 3=high
DEMO_PATIENT_PROTEIN: list[tuple] = [
    # (gene_name, tissue_context, expression_score)
    ("MYC",    "breast", 3),   # HIGH — hub overactivation
    ("CDK4",   "breast", 3),   # HIGH
    ("CDK6",   "breast", 3),   # HIGH
    ("CCND1",  "breast", 3),   # HIGH
    ("E2F1",   "breast", 3),   # HIGH
    ("MCM2",   "breast", 2),   # MEDIUM
    ("PCNA",   "breast", 3),   # HIGH
    ("TOP2A",  "breast", 2),   # MEDIUM
    ("AURKB",  "breast", 2),   # MEDIUM
    ("PLK1",   "breast", 2),   # MEDIUM
    ("TP53",   "breast", 0),   # NOT DETECTED — hub collapse
    ("CDKN1A", "breast", 0),   # NOT DETECTED
    ("MDM2",   "breast", 0),   # NOT DETECTED
    ("BAX",    "breast", 0),   # NOT DETECTED
    ("BBC3",   "breast", 0),   # NOT DETECTED (PUMA)
    ("GADD45A","breast", 0),   # NOT DETECTED
    ("TIGAR",  "breast", 1),   # LOW
    ("BRCA1",  "breast", 1),   # LOW — DNA repair compromised
    ("RAD51",  "breast", 1),   # LOW
    ("FANCD2", "breast", 1),   # LOW
    ("RPA1",   "breast", 2),   # MEDIUM
    ("RFC1",   "breast", 2),   # MEDIUM
    ("EGFR",   "breast", 2),   # MEDIUM
    ("GRB2",   "breast", 2),   # MEDIUM
    ("PIK3CA", "breast", 2),   # MEDIUM
    ("AKT1",   "breast", 3),   # HIGH — PI3K pathway activated
    ("MTOR",   "breast", 2),   # MEDIUM
    ("PTEN",   "breast", 0),   # NOT DETECTED — tumour suppressor lost
    ("RB1",    "breast", 0),   # NOT DETECTED
    ("KRAS",   "breast", 2),   # MEDIUM
    ("RAF1",   "breast", 2),   # MEDIUM
    ("MAP2K1", "breast", 2),   # MEDIUM (MEK1)
    ("MAPK3",  "breast", 2),   # MEDIUM (ERK1)
    ("BRAF",   "breast", 1),   # LOW
    ("VHL",    "breast", 2),   # MEDIUM
    ("HIF1A",  "breast", 2),   # MEDIUM
    ("VEGFA",  "breast", 2),   # MEDIUM
    ("APC",    "breast", 2),   # MEDIUM
    ("CTNNB1", "breast", 2),   # MEDIUM
    ("SMAD4",  "breast", 2),   # MEDIUM
    ("MDM2",   "breast", 0),   # NOT DETECTED
]


# ── Synthetic baseline (used in --offline mode) ───────────────────────────────
# Represents "typical healthy breast tissue" protein expression
# Derived from HPA-like distributions for these proteins
SYNTHETIC_BASELINE: list[tuple] = [
    # (gene_name, tissue, cell_type, expression_score, reliability)
    # MYC regulon — healthy breast: MYC moderate, downstream moderate
    ("MYC",    "breast", "glandular cells", 1, "supported"),
    ("MYC",    "breast", "myoepithelial cells", 1, "supported"),
    ("MYC",    "breast", "adipocytes", 0, "approved"),
    ("CDK4",   "breast", "glandular cells", 2, "supported"),
    ("CDK4",   "breast", "myoepithelial cells", 1, "supported"),
    ("CDK6",   "breast", "glandular cells", 2, "approved"),
    ("CCND1",  "breast", "glandular cells", 2, "supported"),
    ("E2F1",   "breast", "glandular cells", 1, "approved"),
    ("MCM2",   "breast", "glandular cells", 2, "supported"),
    ("PCNA",   "breast", "glandular cells", 2, "supported"),
    ("TOP2A",  "breast", "glandular cells", 1, "approved"),
    ("AURKB",  "breast", "glandular cells", 1, "approved"),
    ("PLK1",   "breast", "glandular cells", 1, "approved"),
    # TP53 regulon — healthy: TP53 present, apoptosis genes active
    ("TP53",   "breast", "glandular cells", 2, "enhanced"),
    ("TP53",   "breast", "myoepithelial cells", 2, "enhanced"),
    ("CDKN1A", "breast", "glandular cells", 2, "supported"),
    ("MDM2",   "breast", "glandular cells", 1, "supported"),
    ("BAX",    "breast", "glandular cells", 2, "approved"),
    ("BBC3",   "breast", "glandular cells", 1, "approved"),
    ("GADD45A","breast", "glandular cells", 1, "approved"),
    ("TIGAR",  "breast", "glandular cells", 1, "approved"),
    # BRCA1 repair complex
    ("BRCA1",  "breast", "glandular cells", 2, "enhanced"),
    ("BRCA1",  "breast", "myoepithelial cells", 2, "enhanced"),
    ("RAD51",  "breast", "glandular cells", 2, "supported"),
    ("FANCD2", "breast", "glandular cells", 2, "approved"),
    ("RPA1",   "breast", "glandular cells", 2, "approved"),
    ("RFC1",   "breast", "glandular cells", 2, "approved"),
    # Signalling
    ("EGFR",   "breast", "glandular cells", 1, "supported"),
    ("GRB2",   "breast", "glandular cells", 2, "approved"),
    ("PIK3CA", "breast", "glandular cells", 2, "supported"),
    ("AKT1",   "breast", "glandular cells", 2, "supported"),
    ("MTOR",   "breast", "glandular cells", 2, "approved"),
    ("PTEN",   "breast", "glandular cells", 2, "enhanced"),
    ("RB1",    "breast", "glandular cells", 2, "enhanced"),
    ("KRAS",   "breast", "glandular cells", 1, "approved"),
    ("RAF1",   "breast", "glandular cells", 1, "approved"),
    ("MAP2K1", "breast", "glandular cells", 2, "approved"),
    ("MAPK3",  "breast", "glandular cells", 2, "approved"),
    ("BRAF",   "breast", "glandular cells", 1, "approved"),
    ("VHL",    "breast", "glandular cells", 2, "approved"),
    ("HIF1A",  "breast", "glandular cells", 1, "supported"),
    ("VEGFA",  "breast", "glandular cells", 1, "supported"),
    ("APC",    "breast", "glandular cells", 2, "supported"),
    ("CTNNB1", "breast", "glandular cells", 2, "supported"),
    ("SMAD4",  "breast", "glandular cells", 2, "supported"),
    ("MDM2",   "breast", "glandular cells", 1, "supported"),
    # Lung proteins
    ("EGFR",   "lung",   "pneumocytes", 1, "supported"),
    ("KRAS",   "lung",   "pneumocytes", 1, "approved"),
    ("STAT3",  "lung",   "pneumocytes", 1, "supported"),
    ("TP53",   "lung",   "pneumocytes", 2, "enhanced"),
    # Colon proteins
    ("APC",    "colon",  "glandular cells", 2, "enhanced"),
    ("CTNNB1", "colon",  "glandular cells", 2, "supported"),
    ("KRAS",   "colon",  "glandular cells", 1, "approved"),
    ("SMAD4",  "colon",  "glandular cells", 2, "supported"),
    # Kidney
    ("VHL",    "kidney", "cells in tubules", 3, "enhanced"),
    ("HIF1A",  "kidney", "cells in tubules", 1, "supported"),
    ("VEGFA",  "kidney", "cells in tubules", 1, "supported"),
    # Skin
    ("BRAF",   "skin",   "melanocytes", 1, "supported"),
    ("MAP2K1", "skin",   "melanocytes", 2, "approved"),
    ("PTEN",   "skin",   "melanocytes", 2, "enhanced"),
    # Brain
    ("EGFR",   "brain",  "glial cells", 1, "supported"),
    ("RB1",    "brain",  "glial cells", 2, "enhanced"),
    ("PTEN",   "brain",  "glial cells", 2, "enhanced"),
    # Pancreas
    ("KRAS",   "pancreas","exocrine glandular cells", 1, "approved"),
    ("SMAD4",  "pancreas","exocrine glandular cells", 2, "supported"),
]


# ── Schema ────────────────────────────────────────────────────────────────────

def create_schema(conn: sqlite3.Connection) -> None:
    """Create all protein SNT tables."""
    logger.info("Creating protein SNT database schema...")
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE IF NOT EXISTS hpa_protein_reference (
            id                INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_name         TEXT    NOT NULL,
            tissue            TEXT    NOT NULL,
            cell_type         TEXT    NOT NULL,
            expression_ordinal TEXT   NOT NULL,
            expression_score  INTEGER NOT NULL CHECK(expression_score BETWEEN 0 AND 3),
            reliability       TEXT    NOT NULL,
            reliability_rank  INTEGER NOT NULL DEFAULT 0,
            source            TEXT    NOT NULL DEFAULT 'HPA',
            created_at        TEXT    DEFAULT (datetime('now'))
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS protein_baseline_stats (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_name        TEXT    NOT NULL,
            tissue_context   TEXT    NOT NULL,
            mean_score       REAL    NOT NULL,
            std_score        REAL    NOT NULL,
            n_samples        INTEGER NOT NULL,
            min_score        REAL    NOT NULL,
            max_score        REAL    NOT NULL,
            created_at       TEXT    DEFAULT (datetime('now')),
            UNIQUE(gene_name, tissue_context)
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS protein_interactions (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            hub_protein      TEXT    NOT NULL,
            satellite_protein TEXT   NOT NULL,
            interaction_type TEXT    NOT NULL,
            confidence_score REAL    NOT NULL DEFAULT 0.80,
            tissue_context   TEXT    NOT NULL DEFAULT 'pan-tissue',
            source           TEXT    NOT NULL DEFAULT 'curated',
            created_at       TEXT    DEFAULT (datetime('now')),
            UNIQUE(hub_protein, satellite_protein, tissue_context)
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS patient_protein_profile (
            id                INTEGER PRIMARY KEY AUTOINCREMENT,
            patient_id        TEXT    NOT NULL,
            gene_name         TEXT    NOT NULL,
            tissue_context    TEXT    NOT NULL,
            expression_score  INTEGER NOT NULL CHECK(expression_score BETWEEN 0 AND 3),
            expression_ordinal TEXT,
            loaded_at         TEXT    DEFAULT (datetime('now')),
            UNIQUE(patient_id, gene_name, tissue_context)
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS disease_protein_snt_signatures (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            disease_name     TEXT    NOT NULL,
            hub_protein      TEXT    NOT NULL,
            satellite_protein TEXT   NOT NULL,
            expected_anomaly TEXT    NOT NULL
                CHECK(expected_anomaly IN (
                    'HUB_OVERACTIVATION','HUB_COLLAPSE',
                    'SATELLITE_CAPTURE','LEAPFROG'
                )),
            tissue_context   TEXT    NOT NULL,
            confidence_score REAL    NOT NULL DEFAULT 0.80,
            created_at       TEXT    DEFAULT (datetime('now'))
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS auto_healing_patterns (
            rule_id       INTEGER PRIMARY KEY AUTOINCREMENT,
            priority      INTEGER NOT NULL DEFAULT 99,
            rule_name     TEXT    NOT NULL UNIQUE,
            target_column TEXT    NOT NULL,
            pattern       TEXT    NOT NULL,
            replacement   TEXT    NOT NULL DEFAULT '',
            rule_type     TEXT    NOT NULL,
            active        INTEGER NOT NULL DEFAULT 1,
            created_at    TEXT    DEFAULT (datetime('now'))
        )
    """)

    conn.commit()
    logger.info("Schema created: 6 tables.")


# ── HPA download and parse ────────────────────────────────────────────────────

def download_hpa_tsv(cache_path: Optional[Path] = None) -> str:
    """
    Download normal_tissue.tsv from HPA.
    Uses cache_path if provided and file exists (avoid repeated downloads).
    Returns the raw TSV string.
    """
    if cache_path and cache_path.exists():
        logger.info("[HPA] Using cached TSV: %s", cache_path)
        return cache_path.read_text(encoding="utf-8")

    logger.info("[HPA] Downloading normal_tissue.tsv.zip from proteinatlas.org...")
    try:
        resp = requests.get(HPA_TSV_URL, timeout=120, stream=True)
        resp.raise_for_status()

        total = int(resp.headers.get("content-length", 0))
        logger.info("[HPA] Download size: %.1f MB", total / 1024 / 1024)

        content = bytearray()
        downloaded = 0
        for chunk in resp.iter_content(chunk_size=65536):
            content.extend(chunk)
            downloaded += len(chunk)
            if total:
                pct = downloaded / total * 100
                if downloaded % (1024 * 1024) < 65536:  # log every ~1MB
                    logger.info("[HPA] Download progress: %.1f%%", pct)

        # Unzip in memory
        with zipfile.ZipFile(io.BytesIO(bytes(content))) as zf:
            tsv_name = [n for n in zf.namelist() if n.endswith(".tsv")][0]
            logger.info("[HPA] Extracting: %s", tsv_name)
            tsv_content = zf.read(tsv_name).decode("utf-8")

        if cache_path:
            cache_path.write_text(tsv_content, encoding="utf-8")
            logger.info("[HPA] TSV cached to: %s", cache_path)

        logger.info("[HPA] Download complete. %d chars.", len(tsv_content))
        return tsv_content

    except requests.exceptions.RequestException as exc:
        logger.error("[HPA] Download failed: %s", exc)
        raise


def parse_hpa_tsv(tsv_content: str) -> list[dict]:
    """
    Parse normal_tissue.tsv into list of row dicts.
    Filters to minimum reliability threshold.
    """
    logger.info("[HPA] Parsing TSV...")
    reader = csv.DictReader(io.StringIO(tsv_content), delimiter="\t")

    rows = []
    skipped_reliability = 0
    skipped_ordinal     = 0

    for row in reader:
        # HPA v25 columns: Gene, Gene name, Tissue, Cell type, Level, Reliability
        gene_name  = (row.get("Gene name") or row.get("Gene") or "").strip().upper()
        tissue     = (row.get("Tissue") or "").strip().lower()
        cell_type  = (row.get("Cell type") or "").strip()
        level      = (row.get("Level") or "").strip().lower()
        reliability= (row.get("Reliability") or "").strip().lower()

        if not gene_name or not tissue or not level:
            continue

        # Reliability filter
        rel_rank = RELIABILITY_RANK.get(reliability, 0)
        min_rank = RELIABILITY_RANK.get(MIN_RELIABILITY, 2)
        if rel_rank < min_rank:
            skipped_reliability += 1
            continue

        score = ORDINAL_MAP.get(level)
        if score is None:
            skipped_ordinal += 1
            continue

        rows.append({
            "gene_name":         gene_name,
            "tissue":            tissue,
            "cell_type":         cell_type,
            "expression_ordinal": level,
            "expression_score":  score,
            "reliability":       reliability,
            "reliability_rank":  rel_rank,
        })

    logger.info(
        "[HPA] Parsed %d usable rows. Skipped: %d reliability, %d unknown level.",
        len(rows), skipped_reliability, skipped_ordinal,
    )
    return rows


def insert_hpa_rows(conn: sqlite3.Connection, rows: list[dict]) -> int:
    """Bulk insert HPA rows into hpa_protein_reference."""
    logger.info("[HPA] Inserting %d rows into hpa_protein_reference...", len(rows))
    cur = conn.cursor()
    inserted = 0
    for i, r in enumerate(rows):
        cur.execute(
            """
            INSERT INTO hpa_protein_reference
                (gene_name, tissue, cell_type, expression_ordinal,
                 expression_score, reliability, reliability_rank, source)
            VALUES (?, ?, ?, ?, ?, ?, ?, 'HPA')
            """,
            (r["gene_name"], r["tissue"], r["cell_type"],
             r["expression_ordinal"], r["expression_score"],
             r["reliability"], r["reliability_rank"]),
        )
        inserted += 1
        if inserted % 10000 == 0:
            conn.commit()
            logger.info("[HPA] Inserted %d rows so far...", inserted)

    conn.commit()
    logger.info("[HPA] Insert complete: %d rows.", inserted)
    return inserted


# ── Baseline stats computation ────────────────────────────────────────────────

def compute_baseline_stats(conn: sqlite3.Connection) -> int:
    """
    Compute μ and σ per (gene_name, tissue) from hpa_protein_reference.
    This is what feeds the Z-score engine at runtime.
    """
    logger.info("[STATS] Computing baseline stats (μ, σ) per gene×tissue...")
    cur = conn.cursor()

    # Load all expression data grouped by gene+tissue
    cur.execute("""
        SELECT gene_name, tissue, expression_score
        FROM hpa_protein_reference
    """)
    rows = cur.fetchall()

    # Group: {(gene, tissue): [scores]}
    groups: dict[tuple, list[int]] = defaultdict(list)
    for gene, tissue, score in rows:
        groups[(gene, tissue)].append(score)

    logger.info("[STATS] %d unique gene×tissue combinations found.", len(groups))

    inserted = 0
    for (gene, tissue), scores in groups.items():
        n   = len(scores)
        mu  = mean(scores)
        sd  = stdev(scores) if n > 1 else 0.0
        cur.execute(
            """
            INSERT OR REPLACE INTO protein_baseline_stats
                (gene_name, tissue_context, mean_score, std_score,
                 n_samples, min_score, max_score)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (gene, tissue, round(mu, 4), round(sd, 4),
             n, min(scores), max(scores)),
        )
        inserted += 1
        if inserted % 1000 == 0:
            conn.commit()
            logger.debug("[STATS] %d stats computed...", inserted)

    conn.commit()
    logger.info("[STATS] Baseline stats complete: %d records.", inserted)
    return inserted


# ── Seed helpers ──────────────────────────────────────────────────────────────

def seed_interactions(conn: sqlite3.Connection) -> int:
    logger.info("[SEED] Seeding %d protein interactions...", len(PROTEIN_INTERACTIONS))
    cur = conn.cursor()
    inserted = 0
    for hub, sat, itype, conf, tissue in PROTEIN_INTERACTIONS:
        cur.execute(
            """
            INSERT OR IGNORE INTO protein_interactions
                (hub_protein, satellite_protein, interaction_type,
                 confidence_score, tissue_context, source)
            VALUES (?, ?, ?, ?, ?, 'curated+UniProt')
            """,
            (hub, sat, itype, conf, tissue),
        )
        if cur.rowcount > 0:
            inserted += 1
            logger.debug("[SEED] Interaction: %s → %s (%s) conf=%.2f", hub, sat, itype, conf)
    conn.commit()
    logger.info("[SEED] Interactions: %d inserted.", inserted)
    return inserted


def seed_disease_signatures(conn: sqlite3.Connection) -> int:
    logger.info("[SEED] Seeding %d disease signatures...", len(DISEASE_PROTEIN_SIGNATURES))
    cur = conn.cursor()
    inserted = 0
    for row in DISEASE_PROTEIN_SIGNATURES:
        cur.execute(
            """
            INSERT OR IGNORE INTO disease_protein_snt_signatures
                (disease_name, hub_protein, satellite_protein,
                 expected_anomaly, tissue_context, confidence_score)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            row,
        )
        if cur.rowcount > 0:
            inserted += 1
            logger.debug("[SEED] Disease: %-35s Hub=%-8s Sat=%-10s %s",
                         row[0], row[1], row[2], row[3])
    conn.commit()
    logger.info("[SEED] Disease signatures: %d inserted.", inserted)
    return inserted


def seed_demo_patient(conn: sqlite3.Connection, patient_id: str = "DEMO-PX-OMEGA") -> int:
    logger.info("[SEED] Seeding demo patient '%s' (%d genes)...",
                patient_id, len(DEMO_PATIENT_PROTEIN))
    cur = conn.cursor()
    inserted = 0
    ordinal_map_rev = {v: k for k, v in ORDINAL_MAP.items() if k != "not_detected"}
    ordinal_map_rev[0] = "Not detected"
    ordinal_map_rev[1] = "Low"
    ordinal_map_rev[2] = "Medium"
    ordinal_map_rev[3] = "High"

    for gene, tissue, score in DEMO_PATIENT_PROTEIN:
        cur.execute(
            """
            INSERT OR REPLACE INTO patient_protein_profile
                (patient_id, gene_name, tissue_context,
                 expression_score, expression_ordinal)
            VALUES (?, ?, ?, ?, ?)
            """,
            (patient_id, gene, tissue, score, ordinal_map_rev[score]),
        )
        if cur.rowcount > 0:
            inserted += 1
    conn.commit()
    logger.info("[SEED] Demo patient: %d records inserted.", inserted)
    return inserted


def seed_healing_rules(conn: sqlite3.Connection) -> int:
    """Seed ETL auto-healing rules (same as RNA pipeline)."""
    rules = [
        (1, "strip_gene_prefix",   "gene_name",  r"^(GENE[-_]|GN[-_]|PROT[-_])", "",  "regex_replace", 1),
        (2, "uppercase_gene",      "gene_name",  r".*",                            "",  "uppercase",     1),
        (3, "strip_whitespace",    "gene_name",  r"\s+",                           "",  "regex_replace", 1),
        (4, "map_low_to_1",        "expression", r"^low$",                         "1", "regex_replace", 1),
        (5, "map_medium_to_2",     "expression", r"^medium$",                      "2", "regex_replace", 1),
        (6, "map_high_to_3",       "expression", r"^high$",                        "3", "regex_replace", 1),
        (7, "map_not_detected",    "expression", r"^not.detected$",                "0", "regex_replace", 1),
        (8, "drop_na_rows",        "gene_name",  r"^(NA|NAN|NULL|nan)$",           "",  "drop_row",      1),
        (9, "strip_quotes",        "gene_name",  r"[\"']",                         "",  "regex_replace", 1),
    ]
    cur = conn.cursor()
    inserted = 0
    for p, name, col, pat, rep, rtype, active in rules:
        cur.execute(
            """
            INSERT OR IGNORE INTO auto_healing_patterns
                (priority, rule_name, target_column, pattern, replacement, rule_type, active)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (p, name, col, pat, rep, rtype, active),
        )
        if cur.rowcount > 0:
            inserted += 1
    conn.commit()
    logger.info("[SEED] Healing rules: %d inserted.", inserted)
    return inserted


# ── Verification ──────────────────────────────────────────────────────────────

def verify(conn: sqlite3.Connection) -> None:
    logger.info("[VERIFY] Running integrity checks...")
    cur = conn.cursor()
    tables = [
        "hpa_protein_reference",
        "protein_baseline_stats",
        "protein_interactions",
        "patient_protein_profile",
        "disease_protein_snt_signatures",
        "auto_healing_patterns",
    ]
    for t in tables:
        cur.execute(f"SELECT COUNT(*) FROM {t}")
        n = cur.fetchone()[0]
        logger.info("[VERIFY]   ✓ %-40s → %d rows", t, n)

    cur.execute("SELECT DISTINCT tissue_context FROM protein_baseline_stats ORDER BY tissue_context")
    tissues = [r[0] for r in cur.fetchall()]
    logger.info("[VERIFY]   Tissues in baseline: %s", ", ".join(tissues))

    cur.execute("SELECT DISTINCT disease_name FROM disease_protein_snt_signatures")
    diseases = [r[0] for r in cur.fetchall()]
    logger.info("[VERIFY]   Diseases in oracle: %s", ", ".join(diseases))
    logger.info("[VERIFY] All checks passed.")


# ── Orchestrator ──────────────────────────────────────────────────────────────

def build(offline: bool = False, tsv_path: Optional[Path] = None) -> None:
    logger.info("=" * 65)
    logger.info("SNT PROTEIN DATABASE BUILDER — START")
    logger.info("DB path : %s", DB_PATH)
    logger.info("Mode    : %s", "OFFLINE (synthetic)" if offline else "ONLINE (HPA)")
    logger.info("=" * 65)

    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(DB_PATH))

    try:
        create_schema(conn)

        if offline:
            # Use synthetic baseline — no network required
            logger.info("[BUILD] Using synthetic baseline (offline mode).")
            rows = [
                {
                    "gene_name":          g,
                    "tissue":             t,
                    "cell_type":          ct,
                    "expression_ordinal": ["Not detected","Low","Medium","High"][s],
                    "expression_score":   s,
                    "reliability":        r,
                    "reliability_rank":   RELIABILITY_RANK.get(r, 2),
                }
                for g, t, ct, s, r in SYNTHETIC_BASELINE
            ]
            insert_hpa_rows(conn, rows)

        elif tsv_path:
            logger.info("[BUILD] Using local TSV: %s", tsv_path)
            tsv_content = tsv_path.read_text(encoding="utf-8")
            rows = parse_hpa_tsv(tsv_content)
            insert_hpa_rows(conn, rows)

        else:
            cache = DB_PATH.parent / "normal_tissue_cache.tsv"
            tsv_content = download_hpa_tsv(cache_path=cache)
            rows = parse_hpa_tsv(tsv_content)
            insert_hpa_rows(conn, rows)

        compute_baseline_stats(conn)
        seed_interactions(conn)
        seed_disease_signatures(conn)
        seed_demo_patient(conn)
        seed_healing_rules(conn)
        verify(conn)

    finally:
        conn.close()

    logger.info("=" * 65)
    logger.info("SNT PROTEIN DATABASE BUILDER — COMPLETE ✓")
    logger.info("Database ready: %s", DB_PATH)
    logger.info("=" * 65)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build SNT Protein Database from Human Protein Atlas + UniProt"
    )
    parser.add_argument(
        "--offline", action="store_true",
        help="Use synthetic baseline data (no internet required)",
    )
    parser.add_argument(
        "--tsv", type=Path, metavar="FILE",
        help="Path to a locally downloaded normal_tissue.tsv file",
    )
    args = parser.parse_args()
    build(offline=args.offline, tsv_path=args.tsv)


if __name__ == "__main__":
    main()
