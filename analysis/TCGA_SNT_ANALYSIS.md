# TCGA SNT Analysis — Corpus Report
**Fractal Core Research | Elán Zainos Corona**  
**Branch:** `genomic-agent-v3`  
**Date:** 2026-06-30  
**Threshold:** |Z| > 2.5 vs cohort baseline

---

## 1. Corpus Summary

| Cohort | Patients | % with Anomaly | % with ≥5 Events | Wall Candidates |
|--------|----------|----------------|------------------|-----------------|
| BRCA   | 1,231    | 29.9%          | 2.1%             | 10 combos       |
| LUAD   | 600      | 31.3%          | 4.8%             | 10 combos       |
| GBM    | 391      | 34.3%          | 4.3%             | 10 combos       |
| COAD   | 524      | 31.5%          | 3.4%             | 10 combos       |
| **TOTAL** | **2,746** | ~31.5%      | ~3.6%            | **40 combos**   |

---

## 2. Top Co-occurrences (Pairs)

### BRCA
Top pairs: `AURKB_UP + PLK1_UP` (15), `BUB1_UP + PLK1_UP` (15), `BRCA2_UP + BUB1_UP` (12)
→ Dominant pattern: mitotic checkpoint genes co-activated

### LUAD
Top pairs: `ATM_UP + BRAF_UP` (11), `ATM_UP + BRCA2_UP` (11), `ATM_UP + PIK3CA_UP` (11)
→ Dominant pattern: ATM hub with DNA repair + oncogenic signaling

### GBM
Top pairs: `BRAF_UP + BRCA2_UP` (8), `PLK1_UP + RAD51_UP` (6), `AURKB_UP + PLK1_UP` (6)
→ Dominant pattern: mitotic instability + DNA repair collapse

### COAD
Top pairs: `APC_UP + BRAF_UP` (9), `APC_UP + KRAS_UP` (9), `APC_UP + PTEN_UP` (9)
→ Dominant pattern: APC (WNT gatekeeper) as central hub co-event

---

## 3. Five-Event Wall Candidates

### BRCA — Top Candidate
`BRCA2_UP | BUB1_UP | FANCD2_UP | PLK1_UP | RAD51_UP` → **4 patients**  
**Interpretation**: Simultaneous overactivation of DNA repair (BRCA2, FANCD2, RAD51) + mitotic checkpoint (BUB1, PLK1) — suggests replication stress response under mitotic override.

### LUAD — Top Candidate (strongest signal)
`ATM_UP | BRAF_UP | BRCA2_UP | PIK3CA_UP | SMAD4_UP` → **9 patients (1.5%)**  
`ATM_UP | BRAF_UP | BRCA2_UP | PTEN_UP | SMAD4_UP` → **9 patients (1.5%)**  
**Interpretation**: ATM (DNA damage sensor) + oncogenic axis (BRAF, PIK3CA) + suppressor loss axis (PTEN) + TGF-β gatekeeper (SMAD4) — multi-pathway simultaneous activation. Strong candidate for pre-invasive signature.

### GBM — Top Candidate
`BRCA1_UP | BUB1_UP | CHEK2_UP | E2F1_UP | TOP2A_UP` → **2 patients**  
**Interpretation**: DNA repair checkpoints (BRCA1, CHEK2) + replication factors (E2F1, TOP2A) + spindle checkpoint (BUB1) — possible replication catastrophe profile.

### COAD — Top Candidate
`APC_UP | ATM_UP | KRAS_UP | PIK3CA_UP | PTEN_UP` → **8 patients (1.5%)**  
**Interpretation**: WNT gatekeeper (APC) + RAS-MAPK (KRAS) + PI3K pathway (PIK3CA, PTEN) + DNA damage sensor (ATM) — classic colorectal oncogenesis plus DNA instability layer.

---

## 4. Cross-Cohort Observations

- **ATM_UP** appears as hub in LUAD (12/29 multi-event) and COAD (8/18 multi-event) → candidate universal sensor
- **BRAF_UP** prominent in GBM, LUAD, COAD → pan-tumor oncogenic driver in multi-event context
- **BRCA2_UP** appears in 3/4 cohorts in multi-event context → possible compensatory overexpression under replication stress
- **PLK1_UP + AURKB_UP** pair: BRCA + GBM → mitotic checkpoint co-activation cross-tumor

---

## 5. Files

| File | Description |
|------|-------------|
| `patient_anomalies.json` | All 2,746 patients with Z-score anomalies |
| `cooccurrence_matrix.json` | Top 100 pairs per cohort |
| `five_event_wall.json` | Candidates with ≥1% frequency filter |
| `five_event_wall_v2.json` | All candidates (no frequency filter) |
| `pipeline_summary.json` | Cohort-level statistics |
| `tcga_snt_results.db` | SQLite: patient_anomalies + wall_candidates |

---

## 6. Methodology

- **Genes analyzed**: 64 hub genes from SNT BASELINE_NETWORK + ACO-A friction genes
- **Anomaly threshold**: |Z| > 2.5 vs cohort-specific baseline (mean ± sd per gene)
- **5-Event Wall**: combinations of ≥5 simultaneous anomalies, empirically derived
- **Data source**: TCGA RNA-seq (STAR counts, tpm_unstranded) via GDC API
- **Pipeline**: Fractal Core Research — genomic-agent-v3 branch
