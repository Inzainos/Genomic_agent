"""
derive_baseline_from_healthy.py — Empirical BASELINE_NETWORK derivation
========================================================================
Derives the healthy-tissue hub-satellite reference ratios used by the SNT
Genomic Topologic Analyzer's Level-1/Level-2 Z-score scanners, replacing the
original hand-calibrated synthetic values in db_builder.BASELINE_NETWORK.

Reference cohort:
  TCGA-BRCA "Solid Tissue Normal" (normal-adjacent) RNA-seq samples,
  STAR-Counts workflow, tpm_unstranded column, via the NIH GDC public API
  (open-access, de-identified). n=40 by default.

For each (hub, satellite) pair already defined in BASELINE_NETWORK, we compute
R = TPM(satellite) / TPM(hub) per healthy sample — the *same* ratio the scanner
computes per patient (agent_logic.run_level1_triage / run_level2_block_scanner)
— then take the cohort mean and population std. std is floored at max(5% of the
mean, 0.01) so the Z-score denominator is always well-defined.

The pairs, gene symbols (including panel aliases like MEK1/4EBP1), and
chromosome assignments are taken verbatim from the existing BASELINE_NETWORK so
only the numeric mean/std change. The NRAS->PI3K pair is left synthetic (PI3K is
a pathway/gene-family label, not a single HGNC symbol).

Reproduce:
    python derive_baseline_from_healthy.py          # uses committed file-id list
    python derive_baseline_from_healthy.py --refresh # re-query GDC for the cohort

Outputs baseline_empirical.json (mean/std/n per pair). The values were then
transcribed into genomic_agent/genomic_database/db_builder.py::BASELINE_NETWORK.
"""
import argparse
import json
import os
import statistics
import sys
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent.resolve()
GENOMIC_AGENT_DIR = HERE.parents[1]
sys.path.insert(0, str(GENOMIC_AGENT_DIR / "genomic_database"))

from db_builder import BASELINE_NETWORK  # noqa: E402

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data/"
IDS_FILE = HERE / "normal_sample_file_ids.txt"
CACHE_DIR = Path(os.environ.get("SNT_TCGA_CACHE", "/tmp/snt_normal_tissue"))

# Panel alias -> HGNC symbol as it appears in the TCGA gene_name column
ALIAS = {
    "MEK1": "MAP2K1", "MEK2": "MAP2K2", "ERK1": "MAPK3", "ERK2": "MAPK1",
    "4EBP1": "EIF4EBP1", "PUMA": "BBC3", "S6K1": "RPS6KB1",
    "TGFb1": "TGFB1", "TGFb2": "TGFB2",
}


def query_normal_file_ids(n: int) -> list[str]:
    filt = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Solid Tissue Normal"]}},
        ],
    }
    params = {"filters": json.dumps(filt), "fields": "file_id", "format": "JSON", "size": str(n)}
    url = GDC_FILES + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as r:
        data = json.load(r)
    return [h["file_id"] for h in data["data"]["hits"]]


def download(fid: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out = CACHE_DIR / f"{fid}.tsv"
    if out.exists() and out.stat().st_size > 100_000:
        return out
    with urllib.request.urlopen(GDC_DATA + fid, timeout=120) as r:
        out.write_bytes(r.read())
    return out


def load_sample(path: Path, wanted: set[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for line in path.read_text().splitlines():
        if not line or line[0] == "#" or line.startswith(("gene_id", "N_")):
            continue
        p = line.split("\t")
        if len(p) >= 7 and p[1] in wanted:
            try:
                out[p[1]] = float(p[6])  # tpm_unstranded
            except ValueError:
                pass
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="re-query GDC for the file-id cohort")
    ap.add_argument("--n", type=int, default=40)
    args = ap.parse_args()

    if args.refresh or not IDS_FILE.exists():
        ids = query_normal_file_ids(args.n)
        IDS_FILE.write_text("\n".join(ids) + "\n")
    else:
        ids = [x.strip() for x in IDS_FILE.read_text().splitlines() if x.strip()]
    print(f"Cohort: {len(ids)} healthy TCGA-BRCA normal-adjacent samples")

    panel = {g for hub, sat, *_ in BASELINE_NETWORK for g in (hub, sat)}
    wanted = {ALIAS.get(g, g) for g in panel}
    samples = [load_sample(download(fid), wanted) for fid in ids]

    def tpm(sample, gene):
        return sample.get(ALIAS.get(gene, gene))

    results = []
    for hub, sat, old_m, old_sd, chrom in BASELINE_NETWORK:
        if "PI3K" in (hub, sat):
            results.append({"hub": hub, "sat": sat, "mean": old_m, "std": old_sd,
                            "n": 0, "chrom": chrom, "empirical": False})
            continue
        ratios = []
        for s in samples:
            h, t = tpm(s, hub), tpm(s, sat)
            if h is None or t is None or h <= 0:
                continue
            ratios.append(t / h)
        if len(ratios) < 10:
            results.append({"hub": hub, "sat": sat, "mean": old_m, "std": old_sd,
                            "n": len(ratios), "chrom": chrom, "empirical": False})
            continue
        mean = statistics.mean(ratios)
        sd = max(statistics.pstdev(ratios) if len(ratios) > 1 else 0.0, 0.05 * mean, 0.01)
        results.append({"hub": hub, "sat": sat, "mean": round(mean, 4), "std": round(sd, 4),
                        "n": len(ratios), "chrom": chrom, "empirical": True})

    out = HERE / "baseline_empirical.json"
    out.write_text(json.dumps(results, indent=2))
    emp = sum(1 for r in results if r["empirical"])
    print(f"Pairs: {len(results)} | empirical: {emp} | retained synthetic: {len(results) - emp}")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
