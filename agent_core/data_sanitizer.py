"""
data_sanitizer.py — SNT ETL Heuristic Data Sanitizer
======================================================
Implements the DataSanitizer class, which:

  1. Loads cleaning rules from the `auto_healing_patterns` table
  2. Applies vectorised Regex corrections on raw RNA-seq DataFrames
  3. Normalises gene IDs, collapses duplicates, validates TPM ranges
  4. Returns a clean DataFrame ready for Z-Score analysis
  5. Logs every repair action for full observability

Design principle: the DB is the source of truth for rules.
Adding a new cleaning pattern requires only a DB insert, not a
code change.

Author  : SNT Genomic Analyzer Team
License : MIT
"""

from __future__ import annotations

import logging
import os
import re
import sqlite3
import sys
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Optional

import pandas as pd

# ── Logging ──────────────────────────────────────────────────────────────────
logger = logging.getLogger("SNT.DataSanitizer")

# ── Config ────────────────────────────────────────────────────────────────────
def _resolve_db() -> Path:
    env = os.getenv("SNT_DB_PATH")
    if env:
        return Path(env)
    for candidate in [Path("/data/snt_genomic.db"),
                      Path(__file__).parent / "snt_genomic.db",
                      Path(__file__).parent.parent / "snt_genomic.db"]:
        if candidate.exists():
            return candidate
    return Path(__file__).parent / "snt_genomic.db"

DB_PATH    = _resolve_db()
TPM_MIN    = 0.0
TPM_MAX    = 1_000_000.0   # biological ceiling for TPM


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class HealingRule:
    """A single auto-healing rule loaded from the database."""
    rule_id:      int
    rule_name:    str
    target_column: str        # 'gene_id' | 'tpm_value' | 'both'
    pattern:      str         # regex pattern to match
    replacement:  str         # replacement string (or special keyword)
    rule_type:    str         # 'regex_replace' | 'drop_row' | 'clip_value' | 'uppercase'
    priority:     int         # lower number = applied first
    active:       bool


@dataclass
class SanitizationReport:
    """Summary of all repairs applied during a sanitization run."""
    input_rows:         int   = 0
    output_rows:        int   = 0
    rules_applied:      int   = 0
    repairs_made:       int   = 0
    rows_dropped:       int   = 0
    duplicates_removed: int   = 0
    out_of_range_clipped: int = 0
    warnings:           list[str] = field(default_factory=list)

    def summary(self) -> str:
        return (
            f"Input={self.input_rows} | Output={self.output_rows} | "
            f"Rules={self.rules_applied} | Repairs={self.repairs_made} | "
            f"Dropped={self.rows_dropped} | Deduped={self.duplicates_removed} | "
            f"Clipped={self.out_of_range_clipped}"
        )


# ── DataSanitizer ─────────────────────────────────────────────────────────────

class DataSanitizer:
    """
    ETL Heuristic Sanitizer for RNA-seq expression data.

    Usage:
        sanitizer = DataSanitizer()
        clean_df, report = sanitizer.sanitize(raw_csv_string)
    """

    def __init__(self) -> None:
        self._rules: list[HealingRule] = []
        self._rules_loaded = False
        logger.info("[SANITIZER] DataSanitizer initialised. DB path: %s", DB_PATH)

    # ── Rule loading ──────────────────────────────────────────────────────────

    def load_rules(self, force_reload: bool = False) -> list[HealingRule]:
        """
        Load active healing rules from the database, ordered by priority.
        Cached after first load unless force_reload=True.
        """
        if self._rules_loaded and not force_reload:
            logger.debug("[SANITIZER] Using cached rules (%d rules).", len(self._rules))
            return self._rules

        logger.info("[SANITIZER] Loading auto_healing_patterns from DB...")
        if not DB_PATH.exists():
            logger.warning("[SANITIZER] DB not found at %s. Running with zero rules.", DB_PATH)
            return []

        conn = sqlite3.connect(str(DB_PATH))
        conn.row_factory = sqlite3.Row
        try:
            rows = conn.execute(
                """
                SELECT rule_id, rule_name, target_column, pattern,
                       replacement, rule_type, priority, active
                FROM auto_healing_patterns
                WHERE active = 1
                ORDER BY priority ASC
                """
            ).fetchall()

            self._rules = [
                HealingRule(
                    rule_id=r["rule_id"],
                    rule_name=r["rule_name"],
                    target_column=r["target_column"],
                    pattern=r["pattern"],
                    replacement=r["replacement"],
                    rule_type=r["rule_type"],
                    priority=r["priority"],
                    active=bool(r["active"]),
                )
                for r in rows
            ]
            self._rules_loaded = True
            logger.info("[SANITIZER] %d active rules loaded.", len(self._rules))
            for rule in self._rules:
                logger.debug(
                    "[SANITIZER]   Rule #%d [P%d] %-30s type=%-15s col=%s",
                    rule.rule_id, rule.priority, rule.rule_name,
                    rule.rule_type, rule.target_column,
                )
            return self._rules
        finally:
            conn.close()

    # ── Core sanitization pipeline ────────────────────────────────────────────

    def sanitize(
        self,
        raw_input: str,
        patient_id: str = "UNKNOWN",
    ) -> tuple[pd.DataFrame, SanitizationReport]:
        """
        Full ETL pipeline for raw RNA-seq CSV input.

        Args:
            raw_input : Raw CSV string (gene_id, tpm_value columns)
            patient_id: Used in log messages for traceability

        Returns:
            (clean_df, report)
              clean_df has columns: gene_id (str), tpm_value (float)
        """
        report = SanitizationReport()
        rules  = self.load_rules()
        report.rules_applied = len(rules)

        logger.info(
            "[SANITIZER] Starting sanitization for patient '%s' (%d chars input).",
            patient_id, len(raw_input),
        )

        # ── Stage 1: Parse raw CSV ────────────────────────────────────────────
        df = self._parse_csv(raw_input, report)
        if df is None or df.empty:
            logger.error("[SANITIZER] CSV parsing failed or empty. Aborting.")
            return pd.DataFrame(columns=["gene_id", "tpm_value"]), report

        report.input_rows = len(df)
        logger.info("[SANITIZER] Stage 1 — Parsed %d rows.", report.input_rows)

        # ── Stage 2: Apply DB-driven healing rules ────────────────────────────
        df, repairs = self._apply_rules(df, rules, report)
        report.repairs_made = repairs
        logger.info("[SANITIZER] Stage 2 — %d repairs applied.", repairs)

        # ── Stage 3: Uppercase normalisation ─────────────────────────────────
        before = df["gene_id"].copy()
        df["gene_id"] = df["gene_id"].str.strip().str.upper()
        changed = (df["gene_id"] != before.str.strip().str.upper()).sum()
        if changed:
            logger.debug("[SANITIZER] Stage 3 — Uppercased %d gene IDs.", changed)

        # ── Stage 4: TPM range validation ─────────────────────────────────────
        df, clipped = self._clip_tpm(df)
        report.out_of_range_clipped = clipped
        logger.info("[SANITIZER] Stage 4 — %d TPM values clipped.", clipped)

        # ── Stage 5: Drop invalid rows ────────────────────────────────────────
        before_drop = len(df)
        df = df[df["gene_id"].str.len() > 0]
        df = df[df["tpm_value"].notna()]
        df = df[df["tpm_value"] >= TPM_MIN]
        dropped = before_drop - len(df)
        report.rows_dropped += dropped
        if dropped:
            logger.warning("[SANITIZER] Stage 5 — Dropped %d invalid rows.", dropped)

        # ── Stage 6: Deduplicate ──────────────────────────────────────────────
        before_dedup = len(df)
        df = df.groupby("gene_id", as_index=False)["tpm_value"].mean()
        deduped = before_dedup - len(df)
        report.duplicates_removed = deduped
        if deduped:
            logger.info(
                "[SANITIZER] Stage 6 — %d duplicate gene IDs collapsed (mean TPM).", deduped
            )

        report.output_rows = len(df)
        logger.info(
            "[SANITIZER] Sanitization complete. %s", report.summary()
        )
        return df, report

    # ── Private helpers ───────────────────────────────────────────────────────

    def _parse_csv(
        self, raw: str, report: SanitizationReport
    ) -> Optional[pd.DataFrame]:
        """Parse raw CSV, auto-detecting delimiter and header."""
        try:
            # Try comma first
            df = pd.read_csv(StringIO(raw), header=None, names=["gene_id", "tpm_value"])

            # Drop header row if first row contains non-numeric tpm
            if not pd.api.types.is_numeric_dtype(df["tpm_value"]):
                try:
                    df["tpm_value"] = pd.to_numeric(df["tpm_value"], errors="coerce")
                except Exception:
                    pass
                # Drop rows where tpm is still NaN (likely the header)
                pre = len(df)
                df = df[df["tpm_value"].notna()]
                dropped = pre - len(df)
                if dropped:
                    logger.info("[SANITIZER] Parse — Dropped %d header/non-numeric rows.", dropped)
                    report.rows_dropped += dropped

            df["tpm_value"] = pd.to_numeric(df["tpm_value"], errors="coerce")
            df["gene_id"]   = df["gene_id"].astype(str)
            logger.debug("[SANITIZER] CSV parsed successfully. Shape: %s", df.shape)
            return df

        except Exception as exc:
            logger.error("[SANITIZER] CSV parse error: %s", exc, exc_info=True)
            report.warnings.append(f"CSV parse error: {exc}")
            return None

    def _apply_rules(
        self,
        df: pd.DataFrame,
        rules: list[HealingRule],
        report: SanitizationReport,
    ) -> tuple[pd.DataFrame, int]:
        """Apply all DB-driven healing rules in priority order."""
        total_repairs = 0

        for rule in rules:
            logger.debug(
                "[SANITIZER] Applying rule #%d '%s' (type=%s, col=%s)",
                rule.rule_id, rule.rule_name, rule.rule_type, rule.target_column,
            )

            try:
                if rule.rule_type == "regex_replace":
                    df, n = self._apply_regex_replace(df, rule)
                    total_repairs += n

                elif rule.rule_type == "drop_row":
                    df, n = self._apply_drop_row(df, rule)
                    report.rows_dropped += n
                    total_repairs += n

                elif rule.rule_type == "uppercase":
                    # Already handled globally in Stage 3
                    pass

                elif rule.rule_type == "clip_value":
                    # Already handled globally in Stage 4
                    pass

                else:
                    logger.warning(
                        "[SANITIZER] Unknown rule type '%s' for rule #%d. Skipping.",
                        rule.rule_type, rule.rule_id,
                    )

            except re.error as exc:
                msg = f"Rule #{rule.rule_id} regex error: {exc}"
                logger.error("[SANITIZER] %s", msg)
                report.warnings.append(msg)

            except Exception as exc:
                msg = f"Rule #{rule.rule_id} unexpected error: {exc}"
                logger.error("[SANITIZER] %s", msg, exc_info=True)
                report.warnings.append(msg)

        return df, total_repairs

    def _apply_regex_replace(
        self, df: pd.DataFrame, rule: HealingRule
    ) -> tuple[pd.DataFrame, int]:
        """Replace regex matches in the target column."""
        col = rule.target_column if rule.target_column != "both" else "gene_id"
        if col not in df.columns:
            logger.warning("[SANITIZER] Column '%s' not in DataFrame. Skipping rule.", col)
            return df, 0

        original = df[col].astype(str).copy()
        df[col] = df[col].astype(str).str.replace(
            rule.pattern, rule.replacement, regex=True
        )
        changed = (df[col] != original).sum()
        if changed:
            logger.info(
                "[SANITIZER]   Rule '%s' → repaired %d cells in col '%s'.",
                rule.rule_name, changed, col,
            )
        return df, int(changed)

    def _apply_drop_row(
        self, df: pd.DataFrame, rule: HealingRule
    ) -> tuple[pd.DataFrame, int]:
        """Drop rows where the target column matches the pattern."""
        col = rule.target_column
        if col not in df.columns:
            return df, 0
        mask   = df[col].astype(str).str.contains(rule.pattern, regex=True, na=False)
        n_drop = mask.sum()
        if n_drop:
            logger.warning(
                "[SANITIZER]   Rule '%s' → dropping %d rows matching pattern '%s'.",
                rule.rule_name, n_drop, rule.pattern,
            )
        df = df[~mask].copy()
        return df, int(n_drop)

    def _clip_tpm(self, df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
        """Clip TPM values to [TPM_MIN, TPM_MAX]. Returns (df, n_clipped)."""
        before = df["tpm_value"].copy()
        df["tpm_value"] = df["tpm_value"].clip(lower=TPM_MIN, upper=TPM_MAX)
        clipped = (df["tpm_value"] != before).sum()
        if clipped:
            logger.warning(
                "[SANITIZER] Clipped %d TPM values outside [%.1f, %.1f].",
                clipped, TPM_MIN, TPM_MAX,
            )
        return df, int(clipped)

    # ── Utility: DataFrame → CSV string ──────────────────────────────────────

    @staticmethod
    def to_csv_string(df: pd.DataFrame) -> str:
        """Serialise clean DataFrame back to CSV for downstream ingestion."""
        return df.to_csv(index=False)
