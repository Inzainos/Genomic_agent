"""
db_builder.py — SNT Genomic Database Builder
============================================
Constructs the three-table SQLite schema that powers the SNT Genomic
Topologic Analyzer Agent.

Tables:
  - baseline_network_reference : healthy tissue hub-satellite ratios (reference)
  - patient_expression          : patient RNA-seq TPM values (runtime ingestion)
  - disease_snt_signatures      : clinical oracle for Level-1 Triage

Run once before launching the agent stack:
    python db_builder.py

Author  : SNT Genomic Analyzer Team
License : MIT
"""

import logging
import os
import sqlite3
import sys
from pathlib import Path

# ── Logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_resolve_log(), mode="w"),
    ],
)
logger = logging.getLogger("SNT.DBBuilder")

# ── Constants ────────────────────────────────────────────────────────────────
def _resolve_db_path() -> Path:
    env = os.getenv("SNT_DB_PATH")
    if env:
        return Path(env)
    # Docker
    docker = Path("/data/snt_genomic.db")
    if Path("/data").exists() and os.access("/data", os.W_OK):
        return docker
    # WSL2 native — put db next to this script
    return Path(__file__).parent / "snt_genomic.db"

DB_PATH = _resolve_db_path()

def _resolve_log() -> str:
    if Path("/data").exists() and os.access("/data", os.W_OK):
        return "/data/db_builder.log"
    return str(Path(__file__).parent / "db_builder.log")

# ── Synthetic reference data ─────────────────────────────────────────────────
# Calibrated to biological realism: MYC regulon + additional hubs
# mean_ratio = TPM(satellite) / TPM(hub) in healthy tissue
# std_dev_ratio = standard deviation across healthy cohort (n=500)

BASELINE_NETWORK = [
    # hub_gene,   satellite_gene,  mean_ratio, std_dev_ratio, chromosome
    ("MYC",       "CDK4",          0.82,        0.09,          "chr8"),
    ("MYC",       "E2F1",          0.75,        0.11,          "chr8"),
    ("MYC",       "CCND1",         0.68,        0.08,          "chr8"),
    ("MYC",       "MCM2",          0.91,        0.07,          "chr8"),
    ("MYC",       "MCM7",          0.88,        0.10,          "chr8"),
    ("MYC",       "PCNA",          0.79,        0.06,          "chr8"),
    ("MYC",       "TOP2A",         0.64,        0.12,          "chr8"),
    ("MYC",       "AURKB",         0.71,        0.09,          "chr8"),
    ("MYC",       "BUB1",          0.55,        0.13,          "chr8"),
    ("MYC",       "PLK1",          0.66,        0.10,          "chr8"),
    # TP53 regulon
    ("TP53",      "CDKN1A",        1.12,        0.14,          "chr17"),
    ("TP53",      "MDM2",          0.43,        0.08,          "chr17"),
    ("TP53",      "BAX",           0.58,        0.09,          "chr17"),
    ("TP53",      "PUMA",          0.49,        0.11,          "chr17"),
    ("TP53",      "GADD45A",       0.37,        0.07,          "chr17"),
    ("TP53",      "TIGAR",         0.31,        0.06,          "chr17"),
    ("TP53",      "DDB2",          0.44,        0.08,          "chr17"),
    ("TP53",      "SESN1",         0.52,        0.10,          "chr17"),
    # BRCA1 regulon
    ("BRCA1",     "RAD51",         0.77,        0.12,          "chr17"),
    ("BRCA1",     "FANCD2",        0.61,        0.09,          "chr17"),
    ("BRCA1",     "RPA1",          0.83,        0.08,          "chr17"),
    ("BRCA1",     "RFC1",          0.70,        0.11,          "chr17"),
    # EGFR regulon
    ("EGFR",      "GRB2",          0.95,        0.07,          "chr7"),
    ("EGFR",      "SOS1",          0.87,        0.09,          "chr7"),
    ("EGFR",      "PIK3CA",        0.74,        0.13,          "chr7"),
    ("EGFR",      "AKT1",          0.81,        0.10,          "chr7"),
    ("EGFR",      "STAT3",         0.66,        0.08,          "chr7"),
    ("EGFR",      "ERK2",          0.90,        0.07,          "chr7"),
    # PIK3CA regulon (chr3)
    ("PIK3CA",    "MTOR",          0.59,        0.11,          "chr3"),
    ("PIK3CA",    "AKT1",          0.88,        0.08,          "chr3"),
    ("PIK3CA",    "S6K1",          0.63,        0.09,          "chr3"),
    ("PIK3CA",    "4EBP1",         0.71,        0.12,          "chr3"),
    # Additional chromosomal coverage for Level-2 block scanner
    ("KRAS",      "RAF1",          0.76,        0.10,          "chr12"),
    ("KRAS",      "MEK1",          0.82,        0.08,          "chr12"),
    ("KRAS",      "ERK1",          0.69,        0.11,          "chr12"),
    ("NRAS",      "RAF1",          0.73,        0.09,          "chr1"),
    ("NRAS",      "PI3K",          0.67,        0.12,          "chr1"),
    ("BRAF",      "MEK1",          0.85,        0.07,          "chr7"),
    ("BRAF",      "MEK2",          0.80,        0.09,          "chr7"),
    ("PTEN",      "AKT1",          0.38,        0.06,          "chr10"),
    ("PTEN",      "MTOR",          0.29,        0.05,          "chr10"),
    ("RB1",       "E2F1",          0.52,        0.08,          "chr13"),
    ("RB1",       "CCND1",         0.45,        0.07,          "chr13"),
    ("VHL",       "HIF1A",         0.34,        0.06,          "chr3"),
    ("VHL",       "VEGFA",         0.41,        0.09,          "chr3"),
    ("CDKN2A",    "CDK4",          0.22,        0.04,          "chr9"),
    ("CDKN2A",    "CDK6",          0.19,        0.03,          "chr9"),
    ("APC",       "CTNNB1",        0.48,        0.08,          "chr5"),
    ("APC",       "TCF4",          0.55,        0.10,          "chr5"),
    ("SMAD4",     "TGFb1",         0.63,        0.11,          "chr18"),
    ("SMAD4",     "TGFb2",         0.58,        0.09,          "chr18"),
]

# ── Disease SNT Signatures (Clinical Oracle) ──────────────────────────────────
# Defines the topological fingerprint of known diseases.
# expected_anomaly: LEAPFROG | SATELLITE_CAPTURE | HUB_COLLAPSE

DISEASE_SNT_SIGNATURES = [
    # Breast Cancer — Basal subtype (TNBC)
    ("Breast_Cancer_Basal_TNBC",    "MYC",    "CDK4",    "LEAPFROG",          "chr8",  0.92),
    ("Breast_Cancer_Basal_TNBC",    "TP53",   "CDKN1A",  "HUB_COLLAPSE",      "chr17", 0.88),
    ("Breast_Cancer_Basal_TNBC",    "BRCA1",  "RAD51",   "SATELLITE_CAPTURE", "chr17", 0.85),
    # Lung Adenocarcinoma
    ("Lung_Adenocarcinoma",         "EGFR",   "STAT3",   "LEAPFROG",          "chr7",  0.90),
    ("Lung_Adenocarcinoma",         "KRAS",   "RAF1",    "SATELLITE_CAPTURE", "chr12", 0.87),
    ("Lung_Adenocarcinoma",         "TP53",   "MDM2",    "HUB_COLLAPSE",      "chr17", 0.82),
    # Colorectal Cancer
    ("Colorectal_Cancer",           "APC",    "CTNNB1",  "HUB_COLLAPSE",      "chr5",  0.89),
    ("Colorectal_Cancer",           "KRAS",   "MEK1",    "SATELLITE_CAPTURE", "chr12", 0.83),
    ("Colorectal_Cancer",           "SMAD4",  "TGFb1",   "HUB_COLLAPSE",      "chr18", 0.78),
    # Melanoma
    ("Melanoma_BRAF_V600E",         "BRAF",   "MEK1",    "SATELLITE_CAPTURE", "chr7",  0.95),
    ("Melanoma_BRAF_V600E",         "BRAF",   "MEK2",    "SATELLITE_CAPTURE", "chr7",  0.91),
    ("Melanoma_BRAF_V600E",         "PTEN",   "AKT1",    "HUB_COLLAPSE",      "chr10", 0.80),
    # Glioblastoma
    ("Glioblastoma_GBM",            "EGFR",   "PIK3CA",  "SATELLITE_CAPTURE", "chr7",  0.88),
    ("Glioblastoma_GBM",            "PTEN",   "MTOR",    "HUB_COLLAPSE",      "chr10", 0.85),
    ("Glioblastoma_GBM",            "RB1",    "CCND1",   "HUB_COLLAPSE",      "chr13", 0.79),
    # Renal Cell Carcinoma
    ("Renal_Cell_Carcinoma",        "VHL",    "HIF1A",   "HUB_COLLAPSE",      "chr3",  0.93),
    ("Renal_Cell_Carcinoma",        "VHL",    "VEGFA",   "HUB_COLLAPSE",      "chr3",  0.90),
    ("Renal_Cell_Carcinoma",        "PTEN",   "AKT1",    "LEAPFROG",          "chr10", 0.75),
    # Pancreatic Ductal Adenocarcinoma
    ("Pancreatic_PDAC",             "KRAS",   "ERK1",    "SATELLITE_CAPTURE", "chr12", 0.91),
    ("Pancreatic_PDAC",             "CDKN2A", "CDK4",    "HUB_COLLAPSE",      "chr9",  0.87),
    ("Pancreatic_PDAC",             "SMAD4",  "TGFb2",   "HUB_COLLAPSE",      "chr18", 0.84),
    # Li-Fraumeni Syndrome (germline TP53)
    ("Li_Fraumeni_Syndrome",        "TP53",   "BAX",     "HUB_COLLAPSE",      "chr17", 0.96),
    ("Li_Fraumeni_Syndrome",        "TP53",   "PUMA",    "HUB_COLLAPSE",      "chr17", 0.94),
    # Hereditary Breast-Ovarian Cancer (BRCA1)
    ("HBOC_BRCA1_Syndrome",         "BRCA1",  "FANCD2",  "HUB_COLLAPSE",      "chr17", 0.91),
    ("HBOC_BRCA1_Syndrome",         "BRCA1",  "RAD51",   "HUB_COLLAPSE",      "chr17", 0.89),
]

# ── Patient Expression — Demo Sample ─────────────────────────────────────────
# Simulates a TNBC patient (MYC overexpression, TP53 loss)
DEMO_PATIENT_EXPRESSION: list[tuple[str, float]] = [
    ("MYC",    450.3),   # Highly overexpressed hub
    ("CDK4",   280.1),   # Elevated satellite → LEAPFROG candidate
    ("E2F1",   310.7),
    ("CCND1",  240.5),
    ("MCM2",   390.2),
    ("MCM7",   375.6),
    ("PCNA",   330.9),
    ("TOP2A",  185.4),
    ("AURKB",  220.3),
    ("BUB1",    95.2),
    ("PLK1",   198.7),
    ("TP53",    22.1),   # Severely suppressed hub → HUB_COLLAPSE
    ("CDKN1A",  18.4),
    ("MDM2",    15.0),
    ("BAX",     10.2),
    ("PUMA",     8.9),
    ("GADD45A",  9.3),
    ("TIGAR",    7.8),
    ("DDB2",    11.2),
    ("SESN1",   12.5),
    ("BRCA1",   30.4),   # Reduced → DNA repair compromise
    ("RAD51",   45.7),
    ("FANCD2",  38.2),
    ("RPA1",    55.3),
    ("RFC1",    49.1),
    ("EGFR",   145.0),
    ("GRB2",   130.2),
    ("SOS1",   120.5),
    ("PIK3CA", 160.3),
    ("AKT1",   175.8),
    ("STAT3",  155.9),
    ("ERK2",   148.4),
    ("PIK3CA", 160.3),
    ("MTOR",   145.7),
    ("KRAS",    88.4),
    ("RAF1",    92.1),
    ("MEK1",    78.3),
    ("ERK1",    85.2),
    ("NRAS",    72.1),
    ("PI3K",    68.9),
    ("BRAF",    95.3),
    ("MEK2",    89.7),
    ("PTEN",    15.2),   # Tumor suppressor silenced
    ("RB1",     19.3),
    ("E2F1",   310.7),   # Freed from RB1 repression
    ("CCND1",  240.5),
    ("VHL",     28.4),
    ("HIF1A",   82.3),
    ("VEGFA",   95.1),
    ("CDKN2A",  12.3),
    ("CDK6",   220.5),
    ("APC",     18.9),
    ("CTNNB1", 215.3),
    ("TCF4",   198.4),
    ("SMAD4",   14.2),
    ("TGFb1",   25.3),
    ("TGFb2",   22.1),
]


# ── Auto-Healing Patterns (ETL Rules) ────────────────────────────────────────
# These rules are loaded at runtime by DataSanitizer from the DB.
# rule_type options: regex_replace | drop_row | uppercase | clip_value
# priority: lower number = applied first

AUTO_HEALING_PATTERNS: list[tuple] = [
    # priority, rule_name, target_column, pattern, replacement, rule_type, active
    (1,  "strip_tpm_prefix",       "gene_id",   r"^TPM[-_]",           "",         "regex_replace", 1),
    (2,  "strip_ens_version",      "gene_id",   r"\.\d+$",             "",         "regex_replace", 1),
    (3,  "replace_comma_decimal",  "tpm_value", r",",                  ".",        "regex_replace", 1),
    (4,  "strip_gene_whitespace",  "gene_id",   r"\s+",                "",         "regex_replace", 1),
    (5,  "remove_special_chars",   "gene_id",   r"[^A-Za-z0-9_\-\.]", "",         "regex_replace", 1),
    (6,  "drop_na_gene_rows",      "gene_id",   r"^(NA|NAN|NULL|nan)$","",         "drop_row",      1),
    (7,  "drop_control_probes",    "gene_id",   r"^CTRL_",             "",         "drop_row",      1),
    (8,  "drop_empty_gene_id",     "gene_id",   r"^\s*$",              "",         "drop_row",      1),
    (9,  "normalise_gene_case",    "gene_id",   r".*",                 "",         "uppercase",     1),
    (10, "clip_tpm_range",         "tpm_value", r".*",                 "",         "clip_value",    1),
    (11, "strip_tpm_suffix",       "gene_id",   r"_TPM$",              "",         "regex_replace", 1),
    (12, "fix_ensg_prefix",        "gene_id",   r"^ensg",              "ENSG",     "regex_replace", 1),
    (13, "drop_deprecated_probes", "gene_id",   r"^DEPRECATED_",       "",         "drop_row",      1),
    (14, "strip_quotes",           "gene_id",   r"[\"']",              "",         "regex_replace", 1),
    (15, "strip_bracket_annot",    "gene_id",   r"\[.*?\]",            "",         "regex_replace", 1),
]


def create_schema(conn: sqlite3.Connection) -> None:
    """Create all four SNT tables if they do not exist."""
    logger.info("Creating SNT database schema...")
    cursor = conn.cursor()

    logger.debug("Creating table: baseline_network_reference")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS baseline_network_reference (
            id             INTEGER PRIMARY KEY AUTOINCREMENT,
            hub_gene       TEXT    NOT NULL,
            satellite_gene TEXT    NOT NULL,
            mean_ratio     REAL    NOT NULL,
            std_dev_ratio  REAL    NOT NULL,
            chromosome     TEXT    NOT NULL DEFAULT 'unknown',
            created_at     TEXT    DEFAULT (datetime('now')),
            UNIQUE(hub_gene, satellite_gene)
        )
    """)

    logger.debug("Creating table: patient_expression")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS patient_expression (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            patient_id  TEXT NOT NULL DEFAULT 'DEMO-PX-001',
            gene_id     TEXT NOT NULL,
            tpm_value   REAL NOT NULL,
            loaded_at   TEXT DEFAULT (datetime('now')),
            UNIQUE(patient_id, gene_id)
        )
    """)

    logger.debug("Creating table: disease_snt_signatures")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS disease_snt_signatures (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            disease_name    TEXT NOT NULL,
            hub_gene_id     TEXT NOT NULL,
            satellite_gene_id TEXT NOT NULL,
            expected_anomaly TEXT NOT NULL
                CHECK(expected_anomaly IN ('LEAPFROG','SATELLITE_CAPTURE','HUB_COLLAPSE')),
            chromosome      TEXT NOT NULL DEFAULT 'unknown',
            confidence_score REAL NOT NULL DEFAULT 0.80,
            created_at      TEXT DEFAULT (datetime('now'))
        )
    """)

    logger.debug("Creating table: auto_healing_patterns")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS auto_healing_patterns (
            rule_id        INTEGER PRIMARY KEY AUTOINCREMENT,
            priority       INTEGER NOT NULL DEFAULT 99,
            rule_name      TEXT    NOT NULL UNIQUE,
            target_column  TEXT    NOT NULL
                CHECK(target_column IN ('gene_id','tpm_value','both')),
            pattern        TEXT    NOT NULL,
            replacement    TEXT    NOT NULL DEFAULT '',
            rule_type      TEXT    NOT NULL
                CHECK(rule_type IN ('regex_replace','drop_row','uppercase','clip_value')),
            active         INTEGER NOT NULL DEFAULT 1,
            created_at     TEXT    DEFAULT (datetime('now'))
        )
    """)

    conn.commit()
    logger.info("Schema created successfully (4 tables).")


def seed_healing_patterns(conn: sqlite3.Connection) -> int:
    """Insert auto-healing ETL rules. Returns rows inserted."""
    logger.info(
        f"Seeding auto_healing_patterns with {len(AUTO_HEALING_PATTERNS)} rules..."
    )
    cursor = conn.cursor()
    inserted = 0
    for row in AUTO_HEALING_PATTERNS:
        priority, rule_name, target_col, pattern, replacement, rule_type, active = row
        try:
            cursor.execute(
                """
                INSERT OR IGNORE INTO auto_healing_patterns
                    (priority, rule_name, target_column, pattern,
                     replacement, rule_type, active)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (priority, rule_name, target_col, pattern, replacement, rule_type, active),
            )
            if cursor.rowcount > 0:
                inserted += 1
                logger.debug(
                    "  [RULE] P%02d %-30s type=%-15s col=%s",
                    priority, rule_name, rule_type, target_col,
                )
        except sqlite3.Error as exc:
            logger.error(f"  [RULE] Insert failed for '{rule_name}': {exc}")

    conn.commit()
    logger.info(f"Healing patterns seeding complete. Inserted {inserted} rules.")
    return inserted


def seed_baseline(conn: sqlite3.Connection) -> int:
    """Insert healthy tissue reference network. Returns rows inserted."""
    logger.info(f"Seeding baseline_network_reference with {len(BASELINE_NETWORK)} records...")
    cursor = conn.cursor()
    inserted = 0
    for row in BASELINE_NETWORK:
        try:
            cursor.execute(
                """
                INSERT OR IGNORE INTO baseline_network_reference
                    (hub_gene, satellite_gene, mean_ratio, std_dev_ratio, chromosome)
                VALUES (?, ?, ?, ?, ?)
                """,
                row,
            )
            if cursor.rowcount > 0:
                inserted += 1
                logger.debug(
                    f"  [BASELINE] Hub={row[0]:<10} Sat={row[1]:<10} "
                    f"μ={row[2]:.3f} σ={row[3]:.3f} chr={row[4]}"
                )
        except sqlite3.Error as exc:
            logger.error(f"  [BASELINE] Insert failed for {row[:2]}: {exc}")

    conn.commit()
    logger.info(f"Baseline seeding complete. Inserted {inserted} records.")
    return inserted


def seed_disease_signatures(conn: sqlite3.Connection) -> int:
    """Insert clinical oracle disease signatures. Returns rows inserted."""
    logger.info(
        f"Seeding disease_snt_signatures with {len(DISEASE_SNT_SIGNATURES)} records..."
    )
    cursor = conn.cursor()
    inserted = 0
    for row in DISEASE_SNT_SIGNATURES:
        try:
            cursor.execute(
                """
                INSERT OR IGNORE INTO disease_snt_signatures
                    (disease_name, hub_gene_id, satellite_gene_id,
                     expected_anomaly, chromosome, confidence_score)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                row,
            )
            if cursor.rowcount > 0:
                inserted += 1
                logger.debug(
                    f"  [DISEASE] {row[0]:<35} Hub={row[1]:<8} "
                    f"Sat={row[2]:<10} Anomaly={row[3]}"
                )
        except sqlite3.Error as exc:
            logger.error(f"  [DISEASE] Insert failed for {row[:3]}: {exc}")

    conn.commit()
    logger.info(f"Disease signature seeding complete. Inserted {inserted} records.")
    return inserted


def seed_demo_patient(conn: sqlite3.Connection, patient_id: str = "DEMO-PX-001") -> int:
    """Insert demo patient RNA-seq data. Returns rows inserted."""
    logger.info(
        f"Seeding demo patient expression for patient_id='{patient_id}' "
        f"({len(DEMO_PATIENT_EXPRESSION)} genes)..."
    )
    cursor = conn.cursor()
    inserted = 0
    for gene_id, tpm in DEMO_PATIENT_EXPRESSION:
        try:
            cursor.execute(
                """
                INSERT OR REPLACE INTO patient_expression
                    (patient_id, gene_id, tpm_value)
                VALUES (?, ?, ?)
                """,
                (patient_id, gene_id, tpm),
            )
            if cursor.rowcount > 0:
                inserted += 1
                logger.debug(f"  [PATIENT] {patient_id} | {gene_id:<12} TPM={tpm:.1f}")
        except sqlite3.Error as exc:
            logger.error(f"  [PATIENT] Insert failed for {gene_id}: {exc}")

    conn.commit()
    logger.info(f"Demo patient seeding complete. Inserted {inserted} records.")
    return inserted


def verify_database(conn: sqlite3.Connection) -> None:
    """Run quick sanity checks and log row counts."""
    logger.info("Running database verification...")
    cursor = conn.cursor()

    tables = ["baseline_network_reference", "patient_expression", "disease_snt_signatures", "auto_healing_patterns"]
    for table in tables:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        count = cursor.fetchone()[0]
        logger.info(f"  ✓ Table [{table}] → {count} rows")

    # Spot-check: list distinct diseases
    cursor.execute("SELECT DISTINCT disease_name FROM disease_snt_signatures ORDER BY disease_name")
    diseases = [r[0] for r in cursor.fetchall()]
    logger.info(f"  ✓ Diseases in oracle: {', '.join(diseases)}")

    # Spot-check: distinct chromosomes in baseline
    cursor.execute(
        "SELECT DISTINCT chromosome FROM baseline_network_reference ORDER BY chromosome"
    )
    chroms = [r[0] for r in cursor.fetchall()]
    logger.info(f"  ✓ Chromosomes in baseline: {', '.join(chroms)}")

    logger.info("Database verification passed.")


def build_database() -> None:
    """Orchestrate the full build pipeline."""
    logger.info("=" * 60)
    logger.info("SNT GENOMIC DATABASE BUILDER — START")
    logger.info(f"Target path: {DB_PATH}")
    logger.info("=" * 60)

    # Ensure data directory exists
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Data directory ensured: {DB_PATH.parent}")

    try:
        conn = sqlite3.connect(str(DB_PATH))
        logger.info(f"Connected to SQLite database: {DB_PATH}")

        create_schema(conn)
        seed_baseline(conn)
        seed_disease_signatures(conn)
        seed_demo_patient(conn)
        seed_healing_patterns(conn)
        verify_database(conn)

        conn.close()
        logger.info("=" * 60)
        logger.info("SNT GENOMIC DATABASE BUILDER — COMPLETE ✓")
        logger.info(f"Database ready at: {DB_PATH}")
        logger.info("=" * 60)

    except sqlite3.Error as exc:
        logger.critical(f"FATAL: Database build failed: {exc}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    build_database()
