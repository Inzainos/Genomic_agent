#!/usr/bin/env python3
"""
Descarga los 263 BRCA faltantes del GDC API y los guarda directamente
como CSV.gz (gene_name + tpm_unstranded) sin guardar el TSV crudo.
"""
import json
import gzip
import csv
import logging
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request
import urllib.error

logging.basicConfig(
    filename='/home/user/workspace/tcga_download_compress.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s'
)
log = logging.getLogger(__name__)

SNT_GENES = {
    'MYC','CDK4','E2F1','CCND1','MCM2','MCM7','PCNA','TOP2A','AURKB','BUB1','PLK1',
    'TP53','CDKN1A','MDM2','BAX','PUMA','GADD45A','TIGAR','DDB2','SESN1',
    'BRCA1','RAD51','FANCD2','RPA1','RFC1',
    'EGFR','GRB2','SOS1','PIK3CA','AKT1','STAT3','ERK2','MTOR','S6K1','4EBP1',
    'KRAS','RAF1','MEK1','ERK1','NRAS','PI3K','BRAF','MEK2','PTEN','RB1',
    'VHL','HIF1A','VEGFA','CDKN2A','CDK6',
    'APC','CTNNB1','TCF4','SMAD4','TGFb1','TGFb2',
    'BRCA2','MLH1','ATM','CHEK2',
}

GDC_BASE = 'https://api.gdc.cancer.gov/data/'
COMP_DIR = Path('/home/user/workspace/tcga_data/compressed/BRCA')
COMP_DIR.mkdir(parents=True, exist_ok=True)

def download_and_compress(entry: dict) -> tuple[str, bool, str]:
    file_id = entry['file_id']
    case_id = entry['case_id']
    out_name = f"{case_id}__{file_id[:8]}.csv.gz"
    out_path = COMP_DIR / out_name

    if out_path.exists():
        return case_id, True, 'already_done'

    url = GDC_BASE + file_id
    retries = 3
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=60) as resp:
                content = resp.read().decode('utf-8', errors='replace')

            rows = []
            for line in content.splitlines():
                line = line.strip()
                if line.startswith('N_') or line.startswith('gene_id') or not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 7:
                    continue
                gene_name = parts[1].strip()
                if gene_name in SNT_GENES:
                    try:
                        tpm = float(parts[6].strip())
                    except:
                        tpm = 0.0
                    rows.append((gene_name, tpm))

            with gzip.open(out_path, 'wt', newline='') as gz:
                w = csv.writer(gz)
                w.writerow(['gene_name', 'tpm_unstranded'])
                w.writerows(rows)

            log.info(f"OK {case_id}: {len(rows)} genes")
            return case_id, True, f'{len(rows)}_genes'

        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 503) and attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            log.error(f"HTTP {e.code} {case_id}: {e}")
            return case_id, False, f'http_{e.code}'
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(1)
                continue
            log.error(f"ERROR {case_id}: {e}")
            return case_id, False, str(e)[:80]

    return case_id, False, 'max_retries'


def main():
    missing = json.loads(open('/home/user/workspace/tcga_data/brca_missing.json').read())
    print(f"Descargando {len(missing)} BRCA faltantes...")
    log.info(f"Iniciando descarga de {len(missing)} BRCA")

    ok_count = 0
    err_count = 0
    errors = []

    with ThreadPoolExecutor(max_workers=6) as ex:
        futures = {ex.submit(download_and_compress, e): e for e in missing}
        for fut in as_completed(futures):
            case_id, ok, reason = fut.result()
            if reason == 'already_done' or ok:
                ok_count += 1
            else:
                err_count += 1
                errors.append({'case_id': case_id, 'reason': reason})

            if (ok_count + err_count) % 50 == 0:
                import subprocess
                r = subprocess.run(['df','-h','/'], capture_output=True, text=True)
                line = r.stdout.split('\n')[1]
                print(f"  {ok_count+err_count}/263 | ok={ok_count} err={err_count} | {line}", flush=True)

    print(f"\n=== DESCARGA COMPLETA: ok={ok_count}, err={err_count} ===")
    log.info(f"FINAL: ok={ok_count}, err={err_count}")

    if errors:
        print(f"Errores: {errors[:10]}")
        with open('/home/user/workspace/tcga_data/brca_download_errors.json', 'w') as f:
            json.dump(errors, f)

    # Verificar total final
    total_comp = sum(len(list((Path('/home/user/workspace/tcga_data/compressed') / c).glob('*.csv.gz')))
                     for c in ['BRCA','LUAD','GBM','COAD'])
    print(f"Total CSV.gz en compressed: {total_comp}")


if __name__ == '__main__':
    main()
