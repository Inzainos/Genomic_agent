"""
app.py — SNT Genomic Topologic Analyzer | Streamlit Interface
=============================================================
Multimodal UI accepting:
  - Text : clinical notes / symptom description
  - File : RNA-seq CSV (gene_id, tpm_value) — runs through DataSanitizer ETL

Full pipeline:
  Input → Guardrails → DataSanitizer ETL → Level 1 Triage →
  Level 2 Deep Scan → LLM Diagnosis → Notifications → Report

Author  : SNT Genomic Analyzer Team
License : MIT
"""

import logging
import os
import sys
import time

import pandas as pd
import streamlit as st

from agent_logic import (
    run_full_analysis,
    verify_guardrails,
    AnalysisReport,
)
from data_sanitizer import DataSanitizer, SanitizationReport

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("/data/streamlit_ui.log", mode="a"),
    ],
)
logger = logging.getLogger("SNT.UI")


# ── Singletons ────────────────────────────────────────────────────────────────
@st.cache_resource
def get_sanitizer() -> DataSanitizer:
    logger.info("[UI] Initialising DataSanitizer singleton...")
    return DataSanitizer()


# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="SNT Genomic Topologic Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
.snt-header {
    background: linear-gradient(135deg, #0d1b2a 0%, #1b263b 50%, #415a77 100%);
    padding: 2rem 2.5rem; border-radius: 12px; margin-bottom: 1.5rem; color: white;
}
.sanity-card {
    background: #1a2a1a; border: 1px solid #2d5a2d; border-radius: 8px;
    padding: 0.8rem 1.2rem; color: #90ee90; font-family: monospace; font-size: 0.82rem;
}
</style>
""", unsafe_allow_html=True)


# ── Session state ─────────────────────────────────────────────────────────────
for key in ("last_report", "last_san_report", "last_clean_df"):
    if key not in st.session_state:
        st.session_state[key] = None


# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## ⚙️ Configuration")

    patient_id = st.text_input(
        "Patient ID",
        value="DEMO-PX-001",
        help="Use DEMO-PX-001 for pre-loaded TNBC synthetic data.",
    )

    st.markdown("---")
    st.markdown("### 🔬 Analysis Thresholds")
    z_threshold = st.slider(
        "Z-Score Anomaly Threshold",
        min_value=1.5, max_value=4.0, value=2.5, step=0.1,
    )
    os.environ["SNT_Z_THRESHOLD"] = str(z_threshold)

    st.markdown("---")
    st.markdown("### 📊 System Status")
    db_path = os.getenv("SNT_DB_PATH", "/data/snt_genomic.db")
    db_ok   = os.path.exists(db_path)
    st.markdown(f"Database:  {'✅ Ready'  if db_ok                            else '❌ Missing'}")
    st.markdown(f"LLM Key:   {'✅ Set'    if os.getenv('OPENROUTER_API_KEY')  else '⚠️  Not set'}")
    st.markdown(f"Telegram:  {'✅ Set'    if os.getenv('TELEGRAM_BOT_TOKEN')  else '⚠️  Not set'}")

    try:
        san   = get_sanitizer()
        rules = san.load_rules()
        st.markdown(f"ETL Rules: ✅ {len(rules)} active")
    except Exception:
        st.markdown("ETL Rules: ⚠️  Error")

    st.markdown("---")
    st.caption("SNT Genomic Analyzer v2.5.0 · Corpus v5 · MIT License")


# ── Header ────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="snt-header">
    <h1 style="margin:0; font-size:2rem;">🧬 SNT Genomic Topologic Analyzer</h1>
    <p style="margin:0.4rem 0 0; opacity:0.8; font-size:1rem;">
        Satellite-Node Topology · Two-Level Scanning · ETL Auto-Healing · SNT v2.5.0
    </p>
</div>
""", unsafe_allow_html=True)

tab_analysis, tab_etl, tab_about, tab_logs = st.tabs(
    ["🔬 Analysis", "🧹 ETL Sanitizer", "ℹ️ About SNT", "📋 Logs"]
)


# ─────────────────────────────────────────────────────────────────────────────
# TAB 1 — ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────
with tab_analysis:
    col_input, col_results = st.columns([1, 1], gap="large")

    with col_input:
        st.markdown("### 📝 Clinical Input")
        clinical_notes = st.text_area(
            "Clinical Notes / Symptoms",
            placeholder=(
                "e.g., 52-year-old female, triple-negative breast mass (2.3 cm), "
                "family history of BRCA1 mutation. FISH shows MYC amplification..."
            ),
            height=160,
        )

        st.markdown("#### 📁 RNA-seq Upload *(optional)*")
        uploaded_file = st.file_uploader(
            "Upload RNA-seq CSV",
            type=["csv"],
            help="Format: gene_id, tpm_value. ETL Auto-Healing runs automatically.",
        )

        if uploaded_file:
            try:
                df_prev = pd.read_csv(uploaded_file, nrows=8, header=None)
                st.dataframe(df_prev, use_container_width=True)
                uploaded_file.seek(0)
            except Exception as exc:
                st.error(f"Preview error: {exc}")

        use_demo = st.checkbox(
            "Use demo patient (DEMO-PX-001)",
            value=True,
            help="Loads pre-seeded TNBC synthetic data. No CSV needed.",
        )
        st.markdown("---")
        run_button = st.button(
            "🚀 Run SNT Analysis", type="primary", use_container_width=True
        )

    with col_results:
        st.markdown("### 📊 Analysis Results")

        if run_button:
            # ── Validation ────────────────────────────────────────────────────
            if not clinical_notes.strip():
                st.error("Please enter clinical notes before running the analysis.")
                st.stop()

            is_safe, reason = verify_guardrails(clinical_notes)
            if not is_safe:
                st.error(f"🛡️ Security Guard: {reason}")
                logger.warning("[UI] Guardrails rejected: %s", reason)
                st.stop()

            active_patient = "DEMO-PX-001" if use_demo else patient_id.strip()
            csv_content    = None

            # ── ETL Sanitization ──────────────────────────────────────────────
            if uploaded_file and not use_demo:
                raw_csv = uploaded_file.read().decode("utf-8")
                logger.info("[UI] Running DataSanitizer on uploaded file...")

                sanitizer = get_sanitizer()
                with st.spinner("🧹 Running ETL Auto-Healing..."):
                    clean_df, san_report = sanitizer.sanitize(raw_csv, active_patient)

                st.session_state.last_san_report = san_report
                st.session_state.last_clean_df   = clean_df

                if clean_df.empty:
                    st.error("ETL produced an empty dataset. Check your CSV format.")
                    st.stop()

                csv_content = DataSanitizer.to_csv_string(clean_df)
                logger.info("[UI] Sanitized CSV: %d genes ready.", len(clean_df))

                if san_report.repairs_made or san_report.rows_dropped:
                    st.warning(
                        f"🧹 ETL repaired **{san_report.repairs_made}** cells, "
                        f"dropped **{san_report.rows_dropped}** rows. "
                        "See **ETL Sanitizer** tab for details."
                    )

            # ── Run pipeline ──────────────────────────────────────────────────
            with st.spinner("Running SNT Two-Level Analysis..."):
                progress = st.progress(0, text="Initialising pipeline...")
                time.sleep(0.2)
                progress.progress(15, text="⚡ Level 1 — Clinical Triage...")

                report: AnalysisReport = run_full_analysis(
                    patient_id=active_patient,
                    clinical_notes=clinical_notes,
                    csv_content=csv_content,
                )

                progress.progress(65, text="🔬 Level 2 — Deep Block Scanner...")
                time.sleep(0.3)
                progress.progress(88, text="🤖 LLM Diagnosis generation...")
                time.sleep(0.3)
                progress.progress(100, text="✅ Complete!")

            st.session_state.last_report = report

            if report.error:
                st.error(f"Analysis error: {report.error}")
                st.stop()

            # ── Metrics ───────────────────────────────────────────────────────
            confirmed = [m for m in report.triage_matches if m.confirmed]
            tentative = [m for m in report.triage_matches if not m.confirmed]
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("✅ Confirmed Diseases", len(confirmed))
            m2.metric("🔎 Tentative Matches",  len(tentative))
            m3.metric("🧬 Orphan Anomalies",   len(report.orphan_anomalies))
            m4.metric("📨 Notifications Sent", len(report.notification_ids))
            st.markdown("---")

            # ── Level 1 ───────────────────────────────────────────────────────
            with st.expander("⚡ Level 1 — Clinical Triage Results", expanded=True):
                if not report.triage_matches:
                    st.info("No triage matches found.")
                else:
                    rows = [
                        {
                            "Disease":    m.disease_name.replace("_", " "),
                            "Hub":        m.hub_gene,
                            "Satellite":  m.satellite_gene,
                            "Chr":        m.chromosome,
                            "Z-Score":    f"{m.detected_z_score:+.3f}",
                            "Type":       m.expected_anomaly,
                            "Status":     "✅ CONFIRMED" if m.confirmed else "⚪ Tentative",
                            "Confidence": f"{m.confidence_score:.0%}",
                        }
                        for m in sorted(
                            report.triage_matches,
                            key=lambda x: abs(x.detected_z_score),
                            reverse=True,
                        )
                    ]
                    st.dataframe(pd.DataFrame(rows), use_container_width=True, height=280)

            # ── Level 2 ───────────────────────────────────────────────────────
            with st.expander("🔬 Level 2 — Orphan Anomalies (Novel Discoveries)", expanded=True):
                if not report.orphan_anomalies:
                    st.success(
                        "No orphan anomalies detected. "
                        "Patient profile aligns with known disease patterns only."
                    )
                else:
                    rows = [
                        {
                            "Chr":       a.chromosome,
                            "Hub":       a.hub_gene,
                            "Satellite": a.satellite_gene,
                            "Z-Score":   f"{a.z_score:+.3f}",
                            "Type":      a.anomaly_type,
                            "R_Patient": f"{a.ratio_patient:.4f}",
                            "R_Healthy": f"{a.ratio_healthy:.4f}",
                        }
                        for a in sorted(
                            report.orphan_anomalies,
                            key=lambda x: abs(x.z_score),
                            reverse=True,
                        )
                    ]
                    st.dataframe(pd.DataFrame(rows), use_container_width=True, height=280)
                    st.info(
                        f"💡 **{len(report.orphan_anomalies)}** orphan anomalies detected — "
                        "patterns NOT in the disease oracle. Possible novel biomarkers."
                    )

            # ── LLM ──────────────────────────────────────────────────────────
            with st.expander("🤖 AI Medical Diagnosis (Claude 3.5 Sonnet)", expanded=True):
                if report.llm_diagnosis.startswith("⚠️"):
                    st.warning(report.llm_diagnosis)
                else:
                    st.markdown(report.llm_diagnosis)

            # ── Notifications ─────────────────────────────────────────────────
            with st.expander("📨 Notification Log"):
                icons = {"jira": "🎫", "slack": "💬", "email": "📧"}
                for svc, nid in report.notification_ids.items():
                    st.markdown(f"{icons.get(svc, '📌')} **{svc.upper()}** — `{nid}`")

        else:
            st.markdown("""
            <div style="text-align:center;padding:3rem;opacity:0.5;">
                <p style="font-size:3rem;">🧬</p>
                <p>Enter clinical notes and click <strong>Run SNT Analysis</strong>.</p>
            </div>
            """, unsafe_allow_html=True)


# ─────────────────────────────────────────────────────────────────────────────
# TAB 2 — ETL SANITIZER
# ─────────────────────────────────────────────────────────────────────────────
with tab_etl:
    st.markdown("### 🧹 ETL Auto-Healing Sanitizer")

    st.markdown("#### Active Healing Rules (live from DB)")
    try:
        san   = get_sanitizer()
        rules = san.load_rules(force_reload=True)
        if rules:
            rule_rows = [
                {
                    "P":         r.priority,
                    "Rule Name": r.rule_name,
                    "Column":    r.target_column,
                    "Type":      r.rule_type,
                    "Pattern":   r.pattern,
                    "Replace":   r.replacement or "(drop/transform)",
                }
                for r in rules
            ]
            st.dataframe(pd.DataFrame(rule_rows), use_container_width=True, height=350)
        else:
            st.warning("No active rules. Has db_builder completed?")
    except Exception as exc:
        st.error(f"Could not load rules: {exc}")

    st.markdown("---")

    # ── Last sanitization report ──────────────────────────────────────────────
    st.markdown("#### Last Sanitization Report")
    san_report: SanitizationReport | None = st.session_state.get("last_san_report")
    clean_df: pd.DataFrame | None         = st.session_state.get("last_clean_df")

    if san_report is None:
        st.info("No sanitization run yet. Upload a CSV in the Analysis tab.")
    else:
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Input Rows",        san_report.input_rows)
        c2.metric("Output Rows",       san_report.output_rows)
        c3.metric("Repairs Made",      san_report.repairs_made)
        c4.metric("Rows Dropped",      san_report.rows_dropped)
        c5.metric("Deduped",           san_report.duplicates_removed)

        if san_report.out_of_range_clipped:
            st.warning(f"⚠️ {san_report.out_of_range_clipped} TPM values clipped.")
        for w in san_report.warnings:
            st.warning(f"⚠️ {w}")

        if clean_df is not None and not clean_df.empty:
            st.dataframe(clean_df.head(20), use_container_width=True)
            st.download_button(
                "⬇️ Download Cleaned CSV",
                data=DataSanitizer.to_csv_string(clean_df),
                file_name="cleaned_expression.csv",
                mime="text/csv",
            )

    st.markdown("---")
    st.markdown("#### 🧪 Live Sanitizer Test Pad")
    test_input = st.text_area(
        "Paste raw CSV to test the ETL pipeline",
        value="TPM-MYC,450.3\nNA,999\nCTRL_PROBE,0.1\nbrca1,30.4\nMYC,450.3\n",
        height=160,
    )
    if st.button("🧹 Run Sanitizer Test"):
        with st.spinner("Sanitizing..."):
            df_test, rep_test = get_sanitizer().sanitize(test_input, "TEST-PAD")
        st.markdown(
            f"<div class='sanity-card'>Input: {rep_test.input_rows} → "
            f"Output: {rep_test.output_rows} | "
            f"Repairs: {rep_test.repairs_made} | "
            f"Dropped: {rep_test.rows_dropped} | "
            f"Deduped: {rep_test.duplicates_removed}</div>",
            unsafe_allow_html=True,
        )
        if not df_test.empty:
            st.dataframe(df_test, use_container_width=True)


# ─────────────────────────────────────────────────────────────────────────────
# TAB 3 — ABOUT
# ─────────────────────────────────────────────────────────────────────────────
with tab_about:
    st.markdown("""
## 🧬 SNT — Shadow Node Theory (Genomic Module)

This application is the genomic implementation of the **Shadow Node Theory (SNT) v2.5.0**,
a universal framework for detecting structural collapse in complex systems through
hub-satellite topological analysis.

### Published Research

- **Preprint (SSRN):** [Shadow Node Theory — SSRN 6418778](https://ssrn.com/abstract=6418778)
- **Repository:** [github.com/Inzainos/The-shadow-Node-Theory](https://github.com/Inzainos/The-shadow-Node-Theory)
- **Zenodo:** [10.5281/zenodo.19446521](https://doi.org/10.5281/zenodo.19446521)
- **Corpus:** v5 · n=721 real cases · 89.3% statistically significant (p=2.5×10⁻⁹⁷)

### SNT Core Hypothesis

> When a master regulator (hub) loses control of its downstream genes (satellites),
> that loss of control **is** the collapse signal — detectable before structural mutations.

| Approach | Traditional Genomics | SNT Topology |
|---|---|---|
| Unit of analysis | Individual mutations | Hub-satellite relationships |
| What it detects | "Wrong letter in code" | "Who stopped controlling whom" |
| Corpus validated | N/A | n=721, ρ=−0.68, p=2.5×10⁻⁹⁷ |
| Novel discovery | Limited to known variants | Orphan anomaly detection |

### Anomaly Types

**🟡 LEAPFROG** — Satellite decouples upward from hub. Z >> +2.5. Oncogene candidate.

**🔵 SATELLITE CAPTURE** — Free gene drawn into hub dependency. Z << -2.5. Tumour suppressor reprogramming.

**🔴 HUB COLLAPSE** — Master regulator loses global control of its regulon. TP53 / BRCA1 / RB1 inactivation.

### Z-Score Formula
```
Z = (R_patient - μ_healthy) / σ_healthy
R = TPM(satellite) / TPM(hub)    |    Threshold: |Z| > 2.5
Corpus: n=721 · n_sig=644 (89.3%) · Spearman ρ=−0.678 · p=2.50×10⁻⁹⁷
```

### ETL Auto-Healing Pipeline
```
Raw CSV → Parse → Regex Rules (DB) → Uppercase → TPM Clip → Drop Invalid → Deduplicate → Clean DF
```
    """)


# ─────────────────────────────────────────────────────────────────────────────
# TAB 4 — LOGS
# ─────────────────────────────────────────────────────────────────────────────
with tab_logs:
    st.markdown("### 📋 Live System Logs")
    log_file = "/data/agent_core.log"
    if os.path.exists(log_file):
        with open(log_file, "r") as f:
            lines = f.readlines()
        last = lines[-150:] if len(lines) > 150 else lines
        level_filter = st.selectbox(
            "Filter", ["ALL", "INFO", "WARNING", "ERROR", "CRITICAL", "DEBUG"], index=0
        )
        if level_filter != "ALL":
            last = [l for l in last if f"| {level_filter}" in l]
        st.text_area("Output", value="".join(last), height=500)
        if st.button("🔄 Refresh"):
            st.rerun()
    else:
        st.info("No log file yet. Run an analysis to populate logs.")
