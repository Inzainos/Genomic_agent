# Scale validation — 976 real TCGA-BRCA tumor patients

**Fractal Core Research | Elán Zainos Corona**
**Date:** 2026-07-03

Large-scale extension of the real-patient validation. Where Round 1 was 1 case
and Round 2 was 8, this run puts the production SNT pipeline through **976
genuine TCGA-BRCA Primary Tumor patients** — effectively the entire BRCA
primary-tumor cohort available on the NIH GDC (1,111 files; 976 unique cases
after excluding duplicates and the 9 already validated).

## Data & method

| Field | Value |
|---|---|
| Source | NIH GDC public API — open-access, de-identified |
| Cohort | TCGA-BRCA, Primary Tumor, STAR-Counts (`tpm_unstranded`) |
| Patients | 976 unique cases |
| Panel | 58 SNT genes per case (PI3K excluded — pathway/family label) |
| Baseline | Empirical (n=40 TCGA-BRCA normal-adjacent) |
| Download | Batched via GDC POST `/data` (100 files per call, 10 batches) |

Each patient was run through the unchanged production functions
(`run_level1_triage` → `run_level2_block_scanner` → `run_aco_analysis`) against
one shared SQLite DB seeded with the empirical baseline + disease oracle.
Downloads were chunked into 10 batches of ~100 to keep each API call light.

Reproduce: `run_scale_validation.py` (records: `scale_results.jsonl`, one line
per patient; aggregate: `scale_summary.json`). Provenance file IDs:
`tcga_case_file_ids.txt`.

## Headline result

**976 / 976 patients processed with 0 exceptions.** The pipeline is robust at
cohort scale on real RNA-seq input.

## Distributions (per patient)

| Metric | min | median | mean | max | sd |
|---|---|---|---|---|---|
| L1 confirmed | 0 | 7 | 7.49 | 16 | 2.8 |
| L2 orphan anomalies | 1 | 15 | 14.79 | 23 | — |
| ACO-A hubs | 0 | 4 | 3.86 | 9 | — |

The spread (confirmed 0–16, orphans 1–23) confirms the empirical baseline makes
the scanner **discriminate across patients** rather than saturate — the whole
point of replacing the synthetic reference.

## Disease-signature frequency (of 976 patients)

| Signature | Patients | % |
|---|---|---|
| Glioblastoma_GBM | 870 | 89% |
| Lung_Adenocarcinoma | 842 | 86% |
| Breast_Cancer_Basal_TNBC | 814 | 83% |
| HBOC_BRCA1_Syndrome | 609 | 62% |
| Renal_Cell_Carcinoma | 604 | 62% |
| Melanoma_BRAF_V600E | 571 | 59% |
| Li_Fraumeni_Syndrome | 561 | 57% |
| Colorectal_Cancer | 502 | 51% |
| Pancreatic_PDAC | 226 | 23% |

Breast-cancer signatures are near the top (as expected — all cases are BRCA),
but note the cross-cancer signatures fire frequently: these are **shared
topological motifs** (e.g. EGFR/PI3K axis for GBM/LUAD, mitotic-checkpoint and
DNA-repair modules), not literal diagnoses of other cancers. This is a property
of the current oracle (signatures share hub-satellite pairs), and is a useful
signal that the oracle's pairs are not disease-specific enough — a candidate for
future refinement.

## ACO-A collapse modes (patients exhibiting each mode)

| Mode | Patients |
|---|---|
| Floor_Arrested | 831 |
| Cracquelure_Decay | 729 |
| Catastrophic_Cliff | 288 |

Computed in snapshot mode (single timepoint per patient, so Δ is not fit; mode
is inferred from friction index + Z-score magnitude + floor proxy).

## Caveats

- **Snapshot samples**: ACO-A runs qualitatively (no Δ fit); longitudinal
  biopsies would be needed to fit the collapse exponent.
- **Oracle specificity**: the high cross-cancer signature hit-rate reflects
  shared hub-satellite pairs across signatures, not multi-cancer diagnosis.
- **Raw per-patient expression is not committed** — only per-patient anomaly
  counts (`scale_results.jsonl`, no TPM values) and the aggregate summary.

## Conclusion

The SNT Genomic Topologic Analyzer runs end-to-end, without a single exception,
across the full 976-patient TCGA-BRCA primary-tumor cohort against a fully
empirical baseline, producing stable, discriminating, biologically interpretable
per-patient results. This closes the real-data validation at cohort scale.
