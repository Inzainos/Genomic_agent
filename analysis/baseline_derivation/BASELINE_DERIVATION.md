# Empirical BASELINE_NETWORK derivation

**Fractal Core Research | Elán Zainos Corona**
**Date:** 2026-07-02

The SNT Genomic Topologic Analyzer scores each patient's hub-satellite
expression ratio R = TPM(satellite)/TPM(hub) against a healthy-tissue reference
(`BASELINE_NETWORK` in `genomic_database/db_builder.py`) via a Z-score. This
directory documents how that reference was upgraded from hand-calibrated
synthetic values to values derived from **real healthy human tissue**.

## Reference cohort

| Field | Value |
|---|---|
| Source | NIH GDC public API (`https://api.gdc.cancer.gov`) — open-access, de-identified |
| Project | TCGA-BRCA |
| Sample type | Solid Tissue Normal (normal-adjacent) |
| Workflow | STAR-Counts, `tpm_unstranded` |
| Cohort size | n = 40 samples |
| Total available | 113 normal-adjacent BRCA samples in GDC |

The 40 GDC file IDs used are listed in `normal_sample_file_ids.txt` for exact
reproducibility.

## Method

For each `(hub, satellite)` pair already defined in `BASELINE_NETWORK`:

1. Compute R = TPM(satellite) / TPM(hub) in each healthy sample — the **same
   ratio the scanner computes per patient** (`agent_logic.run_level1_triage`
   and `run_level2_block_scanner`), skipping samples where the hub is not
   expressed (TPM ≤ 0).
2. `mean_ratio` = cohort mean of R; `std_dev_ratio` = population std of R,
   floored at `max(5% of mean, 0.01)` so the Z-score denominator is always
   well-defined.
3. Panel aliases are mapped to HGNC symbols for lookup in the TCGA `gene_name`
   column (MEK1→MAP2K1, MEK2→MAP2K2, ERK1→MAPK3, ERK2→MAPK1, 4EBP1→EIF4EBP1,
   PUMA→BBC3, S6K1→RPS6KB1, TGFb1→TGFB1, TGFb2→TGFB2).

Pairs, gene symbols, and chromosome assignments are taken verbatim from the
existing `BASELINE_NETWORK`; only the numeric mean/std change.

**Result:** 50/51 pairs derived empirically (all using the full n=40 cohort).
The single `NRAS→PI3K` pair retains its synthetic value — PI3K is a
pathway/gene-family label, not a single HGNC gene symbol, so it cannot be
resolved in the per-gene TSV (its catalytic subunit PIK3CA is captured
separately).

## Files

| File | Description |
|---|---|
| `derive_baseline_from_healthy.py` | Reproducible derivation script |
| `normal_sample_file_ids.txt` | The 40 GDC file IDs of the healthy cohort |
| `baseline_empirical.json` | Derived mean/std/n per pair (transcribed into `db_builder.py`) |

Raw per-sample TSVs (~4 MB each) are **not** committed; the script re-downloads
them from GDC into a local cache on demand.

## Reproduce

```bash
cd genomic_agent/analysis/baseline_derivation
python derive_baseline_from_healthy.py            # uses the committed file-id list
python derive_baseline_from_healthy.py --refresh  # re-query GDC for a fresh cohort
```

The derivation is deterministic given the file-id list: re-running against the
cached cohort reproduces `baseline_empirical.json` exactly (0 mismatches vs the
values baked into `db_builder.py`).

## Effect on the real-patient validation

Re-running `analysis/real_patient_validation/` against the empirical baseline
made the scanner appropriately selective on the real TCGA-BH-A18H tumor sample:

| Metric | Synthetic baseline | Empirical baseline |
|---|---|---|
| Level 1 confirmed | 23 | 10 |
| Level 2 orphan anomalies | 23 | 14 |
| ACO-A hubs processed | 15 | 6 |

Z-scores fell into a biologically plausible range (e.g. `EGFR→STAT3` from
Z=+434.7 to Z=+38); the large values that remain reflect genuine tumor
dysregulation.

## Caveat

The reference is derived from breast-tissue normal-adjacent samples, so it is
most appropriate for breast-cancer patients. Other tissue types would ideally
use a tissue-matched healthy cohort (the same script generalizes by changing
the `project_id`/`sample_type` filter). Normal-adjacent tissue is also not
identical to fully healthy tissue (it can carry field-effect changes), but it
is the standard practical control for TCGA-based analyses.
