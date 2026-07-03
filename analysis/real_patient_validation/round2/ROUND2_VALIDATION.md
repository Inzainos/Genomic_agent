# Round 2 — Batch validation on 8 real TCGA tumor patients

**Fractal Core Research | Elán Zainos Corona**
**Date:** 2026-07-02

Second round of the real-patient validation. Where Round 1 ran a single case
against the (then synthetic) baseline, Round 2 runs a **batch of 8 genuine
TCGA-BRCA Primary Tumor cases** against the **empirical** healthy-tissue baseline
(n=40 normal-adjacent samples) now shipping in `db_builder.BASELINE_NETWORK`.

## Data

| Field | Value |
|---|---|
| Source | NIH GDC public API — open-access, de-identified |
| Cohort | TCGA-BRCA, Primary Tumor, STAR-Counts (`tpm_unstranded`) |
| Patients | 8 (barcodes + GDC file IDs in `tumor_round2_file_ids.txt`) |
| Panel | 58/59 SNT genes per case (PI3K excluded — pathway/family label) |
| Baseline | Empirical (n=40 TCGA-BRCA normal-adjacent) |

Cases: TCGA-AC-A62X, TCGA-AN-A0FN, TCGA-AR-A0TU, TCGA-BH-A18M, TCGA-D8-A1XO,
TCGA-E2-A14U, TCGA-E2-A15T, TCGA-E9-A1R3.

## Method

For each patient, the production `agent_logic` functions were run unchanged:
`run_level1_triage` → `run_level2_block_scanner` → `run_aco_analysis`, against a
single SQLite DB seeded with the empirical baseline + disease oracle. LLM
diagnosis and notifications were skipped (out of scope for scanner validation).

Reproduce:
```bash
cd genomic_agent/analysis/real_patient_validation/round2
python run_round2_batch.py
```

## Results

| Case | Genes | L1 confirmed | L2 orphans | ACO-A hubs |
|---|---|---|---|---|
| TCGA-AC-A62X | 58 | 8 | 15 | 4 |
| TCGA-AN-A0FN | 58 | 7 | 14 | 3 |
| TCGA-AR-A0TU | 58 | 5 | 18 | 2 |
| TCGA-BH-A18M | 58 | 2 | 12 | 0 |
| TCGA-D8-A1XO | 58 | 4 | 14 | 1 |
| TCGA-E2-A14U | 58 | 4 | 7 | 2 |
| TCGA-E2-A15T | 58 | 7 | 18 | 4 |
| TCGA-E9-A1R3 | 58 | 7 | 12 | 4 |
| **mean** | 58 | **5.5** | **13.75** | **2.5** |

Full per-patient detail (confirmed diseases + ACO-A collapse modes):
`round2_summary.json`. **Exceptions across all 8 patients: 0.**

## Observations

- **Biologically sensible.** All 8 are breast tumors; the breast-cancer
  signatures (`Breast_Cancer_Basal_TNBC`, `HBOC_BRCA1_Syndrome`) recur most
  often, with cross-cancer signatures firing where shared pathways
  (EGFR/PI3K, TP53, mitotic checkpoint) are dysregulated.
- **Per-patient heterogeneity is real, not noise.** Confirmed matches range 2–8
  and orphan anomalies 7–18 across cases — the scanner discriminates between
  tumors rather than flagging everything uniformly, which is the behavior the
  empirical baseline was meant to produce (Round 1 under the synthetic baseline
  saturated at 23 confirmed for its single case).
- **ACO-A modes vary by patient**: mostly `Floor_Arrested` and
  `Catastrophic_Cliff`, with `Cracquelure_Decay` in two cases — computed from
  snapshot data (no per-patient timeseries, so Δ is not fit; mode is inferred
  from friction index + Z-score magnitude + floor proxy).

## Caveat

These are snapshot (single-timepoint) tumor samples, so ACO-A runs in its
qualitative/snapshot mode (no Δ fit). Confirming collapse-mode classification
quantitatively would need longitudinal biopsies. The disease oracle and
baseline are both breast-tissue-appropriate here (all cases are BRCA).

## Conclusion

The SNT pipeline runs end-to-end without error on a batch of 8 genuine TCGA
tumor patients against a fully empirical baseline, producing per-patient,
biologically plausible, and appropriately heterogeneous results — extending the
single-case Round 1 validation to a small cohort.
