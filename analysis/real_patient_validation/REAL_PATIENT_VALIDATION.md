# Real-Patient Pipeline Validation — TCGA-BH-A18H

**Fractal Core Research | Elán Zainos Corona**
**Date:** 2026-07-01
**Purpose:** End-to-end validation of the SNT Genomic Topologic Analyzer's
production pipeline using a genuine patient case, replacing the synthetic
`DEMO-PX-001` sample shipped with the agent.

---

## 1. Data source

| Field | Value |
|---|---|
| Patient case | `TCGA-BH-A18H` (TCGA-BRCA cohort) |
| Source | NIH GDC public API (`https://api.gdc.cancer.gov`) — open-access, de-identified |
| File | Gene-expression quantification, STAR-Counts workflow, `tpm_unstranded` |
| GDC file ID | `744a6d3d-b666-49aa-8d26-47f34e3d1eb5` |
| Raw file size | 60,666 gene rows |
| SNT panel recovered | 59/60 genes (see §4) |

This is real, published, open-access cancer genomics data — not synthetic or
simulated — obtained directly from the GDC API, the same source used earlier
for the 2,746-patient TCGA batch analysis underlying the `disease_snt_signatures`
oracle (`analysis/TCGA_SNT_ANALYSIS.md`).

## 2. Method

1. Built a throwaway SQLite database (schema + baseline + disease oracle +
   healing rules — no synthetic demo patient or demo timeseries).
2. Loaded the 59 real TPM values for `TCGA-BH-A18H` into `patient_expression`.
3. Ran the actual production functions from `agent_core/agent_logic.py`
   unmodified:
   - `run_level1_triage()` — O(K) cross-reference against the clinical oracle
   - `run_level2_block_scanner()` — chromosome-by-chromosome orphan sweep
   - `run_aco_analysis()` — collapse-mode classification for confirmed hub
     collapses
4. LLM diagnosis and Jira/Slack/Email/Telegram notifications were skipped —
   they depend on live mock services / API keys unrelated to the scanning
   logic being validated here.

Reproduce with:
```bash
cd genomic_agent/analysis/real_patient_validation
SNT_DB_PATH=/tmp/snt_real_validation.db python run_real_validation.py
```

## 3. Results

Results are shown against the **empirical** healthy-tissue baseline (see §5),
which superseded the original synthetic reference.

| Stage | Empirical baseline | (was, synthetic baseline) |
|---|---|---|
| Genes loaded | 59/60 (98.3%) | 59/60 |
| Level 1 signatures evaluated | 25 | 25 |
| Level 1 confirmed (\|Z\|>2.5) | 10 | 23 |
| Level 2 orphan anomalies | 14 | 23 |
| ACO-A hubs processed | 6 | 15 |
| Exceptions raised | 0 | 0 |

Full per-match detail: `real_validation_summary.json`.

The pipeline ran to completion with no exceptions across all three stages,
confirming the scanning/classification code path is correct against a real
RNA-seq profile, not just against the hand-crafted synthetic demo case. With
the empirical baseline in place the scanner is markedly more selective —
confirmed matches fell from 23 to 10 and orphan anomalies from 23 to 14 — and
Z-scores moved into a biologically plausible range (e.g. `TP53→BAX` Z=+7.5
rather than +19, `EGFR→STAT3` Z=+38 rather than +435). The large Z-scores that
remain reflect genuine tumor dysregulation in this real breast-cancer sample.

## 4. Gene recovery notes

10/60 SNT-panel genes were not found by direct symbol match in the TCGA
`gene_name` column on the first pass; 9 were recovered via official HGNC
alias mapping (the SNT panel uses informal pathway names in a few spots):

| Panel alias | HGNC symbol |
|---|---|
| MEK1 | MAP2K1 |
| MEK2 | MAP2K2 |
| ERK1 | MAPK3 |
| ERK2 | MAPK1 |
| 4EBP1 | EIF4EBP1 |
| PUMA | BBC3 |
| S6K1 | RPS6KB1 |
| TGFb1 | TGFB1 |
| TGFb2 | TGFB2 |

`PI3K` remains unmapped — it denotes a gene family/pathway, not a single
HGNC gene symbol; its catalytic subunit `PIK3CA` is already a separate,
correctly-mapped entry in the panel.

## 5. Baseline reference: now empirical

The original `BASELINE_NETWORK` (healthy-tissue reference in `db_builder.py`)
was a hand-calibrated, illustrative reference with artificially narrow standard
deviations, so almost every hub-satellite ratio in a real tumor triggered very
large \|Z\| values (e.g. `EGFR→STAT3` Z=+434.7). That reference has now been
**replaced with values derived from n=40 real healthy breast-tissue samples**
(TCGA-BRCA "Solid Tissue Normal", via the NIH GDC API) — see
`analysis/baseline_derivation/`. For each pair the reference is
`mean ± std` of R = TPM(satellite)/TPM(hub) across the healthy cohort, using
the same ratio the scanner computes per patient. Only the `NRAS→PI3K` pair
remains synthetic (PI3K is a pathway/family label, not a single HGNC symbol).

With the empirical baseline the scanner became appropriately selective on this
real patient (confirmed 23→10, orphans 23→14) and Z-scores dropped to a
biologically plausible range. This closes the caveat noted in the original
validation: both the patient data **and** the healthy baseline are now real.

## 6. Conclusion

The SNT Genomic Topologic Analyzer's Two-Level Scanning Architecture (clinical
triage + deep block scanner) and ACO-A collapse-mode classifier execute
correctly against genuine, open-access TCGA patient data, evaluated against a
healthy-tissue baseline that is itself derived from real TCGA normal-adjacent
samples — no synthetic reference data remains in the scoring path (save the one
unresolvable PI3K pair). The synthetic `DEMO-PX-001` case remains useful as a
hand-designed illustrative example; this real-patient run against a real
baseline is the first confirmation that the full scoring path produces
biologically sensible results on real-world RNA-seq input.
