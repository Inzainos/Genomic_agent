#!/usr/bin/env python3
"""
Pipeline SNT completo sobre corpus TCGA.
1. Carga todos los CSV.gz por cohorte
2. Calcula Z-scores por gen vs baseline de la cohorte
3. Detecta anomalías topológicas (|Z| > 2.5) en genes hub
4. Genera matriz de co-ocurrencias
5. Identifica candidatos al 5-Event Wall
"""
import gzip
import csv
import json
import logging
import sqlite3
from pathlib import Path
from collections import defaultdict
import statistics
import itertools

logging.basicConfig(
    filename='/home/user/workspace/snt_pipeline.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s'
)
log = logging.getLogger(__name__)

# ── CONFIGURACIÓN ────────────────────────────────────────────────────────────
COMP_DIR = Path('/home/user/workspace/tcga_data/compressed')
OUT_DIR  = Path('/home/user/workspace/tcga_results')
OUT_DIR.mkdir(exist_ok=True)
COHORTS  = ['BRCA', 'LUAD', 'GBM', 'COAD']
Z_THRESH = 2.5   # umbral de anomalía topológica

# Genes hub principales por módulo funcional SNT
HUB_GENES = [
    'MYC','CDK4','E2F1','CCND1',            # módulo proliferación
    'TP53','MDM2','CDKN1A','RB1',           # módulo supresión tumoral
    'BRCA1','RAD51','FANCD2','ATM',         # módulo reparación DNA
    'EGFR','KRAS','BRAF','PTEN',            # módulo señalización
    'AKT1','MTOR','PIK3CA','STAT3',         # módulo mTOR/PI3K
    'VHL','HIF1A','VEGFA',                  # módulo hipoxia
    'APC','CTNNB1','SMAD4','TGFb1',         # módulo WNT/TGF
    'AURKB','BUB1','PLK1','TOP2A',          # módulo mitosis
    'TP53','BRCA2','MLH1','CHEK2',          # genes fricción ACO-A adicionales
]
HUB_GENES = list(dict.fromkeys(HUB_GENES))  # deduplicar preservando orden

FRICTION_GENES = ['TP53','BRCA1','BRCA2','MLH1','ATM','CHEK2','RAD51','FANCD2','RB1','PTEN']


# ── CARGA DE DATOS ────────────────────────────────────────────────────────────
def load_cohort(cohort: str) -> list[dict]:
    """Carga todos los CSV.gz de una cohorte → lista de {case_id, genes: {gene: tpm}}"""
    cohort_dir = COMP_DIR / cohort
    patients = []
    for gz_file in cohort_dir.glob('*.csv.gz'):
        case_id = gz_file.stem.split('__')[0]  # TCGA-XX-XXXX
        genes = {}
        try:
            with gzip.open(gz_file, 'rt') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    gene = row.get('gene_name','').strip()
                    try:
                        tpm = float(row.get('tpm_unstranded', 0))
                    except:
                        tpm = 0.0
                    if gene:
                        genes[gene] = tpm
        except Exception as e:
            log.warning(f"Error leyendo {gz_file.name}: {e}")
            continue
        patients.append({'case_id': case_id, 'cohort': cohort, 'genes': genes})
    log.info(f"[{cohort}] Cargados {len(patients)} pacientes")
    return patients


# ── ESTADÍSTICAS BASELINE ────────────────────────────────────────────────────
def compute_baseline(patients: list[dict]) -> dict:
    """Calcula media y sd por gen sobre todos los pacientes de la cohorte."""
    gene_values = defaultdict(list)
    for p in patients:
        for gene, tpm in p['genes'].items():
            gene_values[gene].append(tpm)

    baseline = {}
    for gene, vals in gene_values.items():
        if len(vals) < 5:
            continue
        mean = statistics.mean(vals)
        try:
            sd = statistics.stdev(vals)
        except:
            sd = 0.0
        baseline[gene] = {'mean': mean, 'sd': sd}
    return baseline


# ── Z-SCORES Y ANOMALÍAS ─────────────────────────────────────────────────────
def detect_anomalies(patient: dict, baseline: dict) -> list[str]:
    """
    Retorna lista de genes hub con |Z| > Z_THRESH.
    Considera anomalía: sobreexpresión (Z > +2.5) o silenciamiento (Z < -2.5).
    """
    anomalies = []
    for gene in HUB_GENES:
        tpm = patient['genes'].get(gene)
        if tpm is None:
            continue
        stats = baseline.get(gene)
        if not stats or stats['sd'] < 1e-6:
            continue
        z = (tpm - stats['mean']) / stats['sd']
        if abs(z) >= Z_THRESH:
            direction = 'UP' if z > 0 else 'DN'
            anomalies.append(f"{gene}_{direction}")
    return anomalies


def compute_friction(patient: dict) -> float:
    """Índice de fricción ACO-A: promedio TPM de genes guardianes / 100"""
    vals = [patient['genes'].get(g, 0.0) for g in FRICTION_GENES]
    vals = [v for v in vals if v > 0]
    if not vals:
        return 0.0
    return statistics.mean(vals) / 100.0


# ── PIPELINE PRINCIPAL ────────────────────────────────────────────────────────
def run_pipeline():
    all_patients = []
    cohort_baselines = {}

    print("Cargando corpus TCGA...", flush=True)
    for cohort in COHORTS:
        patients = load_cohort(cohort)
        baseline = compute_baseline(patients)
        cohort_baselines[cohort] = baseline
        all_patients.extend(patients)
        print(f"  [{cohort}] {len(patients)} pacientes cargados, baseline calculado", flush=True)
        log.info(f"[{cohort}] baseline genes: {len(baseline)}")

    print(f"\nTotal pacientes: {len(all_patients)}", flush=True)

    # ── Detectar anomalías por paciente ─────────────────────────────────────
    print("\nDetectando anomalías topológicas...", flush=True)
    patient_results = []
    for p in all_patients:
        baseline = cohort_baselines[p['cohort']]
        anomalies = detect_anomalies(p, baseline)
        friction = compute_friction(p)
        patient_results.append({
            'case_id': _pseudo_id(p['case_id']),  # SECURITY: pseudonymized barcode
            'cohort': p['cohort'],
            'n_anomalies': len(anomalies),
            'anomalies': anomalies,
            'friction': round(friction, 4),
        })

    # Estadísticas generales
    for cohort in COHORTS:
        coh_pats = [r for r in patient_results if r['cohort'] == cohort]
        mean_anom = statistics.mean([r['n_anomalies'] for r in coh_pats]) if coh_pats else 0
        print(f"  [{cohort}] media anomalías/paciente: {mean_anom:.2f}", flush=True)
        log.info(f"[{cohort}] mean_anomalies={mean_anom:.2f}")

    # ── Guardar resultados por paciente ─────────────────────────────────────
    patient_out = OUT_DIR / 'patient_anomalies.json'
    with open(patient_out, 'w') as f:
        json.dump(patient_results, f, indent=2)
    log.info(f"Guardado patient_anomalies.json ({len(patient_results)} registros)")
    print(f"\nGuardado: {patient_out}", flush=True)

    return patient_results, cohort_baselines


# ── MATRIZ DE CO-OCURRENCIAS ──────────────────────────────────────────────────
def build_cooccurrence_matrix(patient_results: list[dict]) -> dict:
    """
    Cuenta cuántas veces aparece cada par de anomalías juntas en el mismo paciente,
    por cohorte. Solo cuenta pacientes con ≥2 anomalías.
    """
    print("\nConstruyendo matriz de co-ocurrencias...", flush=True)
    cohort_cooc = {}

    for cohort in COHORTS:
        coh_pats = [r for r in patient_results if r['cohort'] == cohort]
        pair_counts = defaultdict(int)
        single_counts = defaultdict(int)
        n_patients = len(coh_pats)

        for p in coh_pats:
            anom = sorted(p['anomalies'])
            for a in anom:
                single_counts[a] += 1
            if len(anom) >= 2:
                for a, b in itertools.combinations(anom, 2):
                    pair_counts[(a, b)] += 1

        # Convertir a lista ordenada por frecuencia
        pairs_sorted = sorted(pair_counts.items(), key=lambda x: x[1], reverse=True)
        cohort_cooc[cohort] = {
            'n_patients': n_patients,
            'single_counts': dict(single_counts),
            'top_pairs': [(list(k), v) for k, v in pairs_sorted[:100]],
        }
        print(f"  [{cohort}] {len(pair_counts)} pares únicos | top: {pairs_sorted[:3] if pairs_sorted else '—'}", flush=True)
        log.info(f"[{cohort}] pairs={len(pair_counts)}")

    cooc_out = OUT_DIR / 'cooccurrence_matrix.json'
    with open(cooc_out, 'w') as f:
        json.dump(cohort_cooc, f, indent=2)
    print(f"Guardado: {cooc_out}", flush=True)
    return cohort_cooc


# ── 5-EVENT WALL ─────────────────────────────────────────────────────────────
def find_five_event_wall(patient_results: list[dict]) -> dict:
    """
    Busca combinaciones de exactamente 5 anomalías que aparezcan con frecuencia
    estadística por cohorte. Solo combina genes con frecuencia individual >5%.
    """
    print("\nBuscando candidatos al 5-Event Wall...", flush=True)
    wall_candidates = {}

    for cohort in COHORTS:
        coh_pats = [r for r in patient_results if r['cohort'] == cohort]
        n = len(coh_pats)
        if n == 0:
            continue

        # Frecuencia individual de cada anomalía
        freq = defaultdict(int)
        for p in coh_pats:
            for a in p['anomalies']:
                freq[a] += 1

        # Solo genes con freq > 5% del cohort
        min_freq = max(3, int(n * 0.03))  # 3% mínimo
        common_anom = sorted([a for a, c in freq.items() if c >= min_freq],
                             key=lambda x: freq[x], reverse=True)
        print(f"  [{cohort}] anomalías comunes (≥{min_freq} pacientes): {len(common_anom)}", flush=True)
        log.info(f"[{cohort}] common_anomalies={len(common_anom)}")

        if len(common_anom) < 5:
            wall_candidates[cohort] = []
            continue

        # Contar co-ocurrencias de combinaciones de 5
        combo5_counts = defaultdict(int)
        for p in coh_pats:
            # Solo anomalías comunes en este paciente
            p_common = [a for a in p['anomalies'] if a in set(common_anom)]
            if len(p_common) >= 5:
                # Limitar a top-20 más frecuentes para reducir explosión combinatoria
                p_top = sorted(p_common, key=lambda x: freq[x], reverse=True)[:20]
                for combo in itertools.combinations(p_top, 5):
                    combo5_counts[combo] += 1

        # Filtrar: aparece en ≥1% de pacientes del cohort
        min_combo = max(2, int(n * 0.01))
        significant = [(list(k), v, round(v/n*100, 2))
                       for k, v in combo5_counts.items() if v >= min_combo]
        significant.sort(key=lambda x: x[1], reverse=True)

        print(f"  [{cohort}] combos 5-eventos significativos: {len(significant)}", flush=True)
        wall_candidates[cohort] = significant[:50]  # top 50 por cohorte
        log.info(f"[{cohort}] wall_candidates={len(significant)}")

        # Mostrar top 5
        for combo, count, pct in significant[:5]:
            print(f"    {' | '.join(combo)} → {count} pacientes ({pct}%)", flush=True)

    wall_out = OUT_DIR / 'five_event_wall.json'
    with open(wall_out, 'w') as f:
        json.dump(wall_candidates, f, indent=2)
    print(f"\nGuardado: {wall_out}", flush=True)
    return wall_candidates


# ── RESUMEN FINAL ─────────────────────────────────────────────────────────────
def generate_summary(patient_results, cooc, wall):
    """Genera resumen estadístico para el reporte."""
    summary = {}
    for cohort in COHORTS:
        coh_pats = [r for r in patient_results if r['cohort'] == cohort]
        n = len(coh_pats)
        if n == 0:
            continue
        n_with_anom = sum(1 for p in coh_pats if len(p['anomalies']) >= 1)
        n_multi = sum(1 for p in coh_pats if len(p['anomalies']) >= 5)
        mean_f = statistics.mean([p['friction'] for p in coh_pats]) if coh_pats else 0
        summary[cohort] = {
            'n_patients': n,
            'n_with_any_anomaly': n_with_anom,
            'pct_with_anomaly': round(n_with_anom/n*100, 1) if n else 0,
            'n_multi_event_5plus': n_multi,
            'pct_multi_event': round(n_multi/n*100, 1) if n else 0,
            'mean_friction': round(mean_f, 4),
            'n_wall_candidates': len(wall.get(cohort, [])),
        }

    summary_out = OUT_DIR / 'pipeline_summary.json'
    with open(summary_out, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nGuardado: {summary_out}", flush=True)
    return summary


# ── GUARDAR EN SQLITE ──────────────────────────────────────────────────────────
def save_to_sqlite(patient_results, wall):
    db_path = OUT_DIR / 'tcga_snt_results.db'
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    c.execute('''CREATE TABLE IF NOT EXISTS patient_anomalies (
        case_id TEXT, cohort TEXT, n_anomalies INTEGER,
        anomalies_json TEXT, friction REAL
    )''')
    c.execute('DELETE FROM patient_anomalies')
    c.executemany(
        'INSERT INTO patient_anomalies VALUES (?,?,?,?,?)',
        [(_pseudo_id(r['case_id']), r['cohort'], r['n_anomalies'],  # SECURITY: pseudonymized
          json.dumps(r['anomalies']), r['friction'])
         for r in patient_results]
    )

    c.execute('''CREATE TABLE IF NOT EXISTS wall_candidates (
        cohort TEXT, events_json TEXT, n_patients INTEGER, pct REAL
    )''')
    c.execute('DELETE FROM wall_candidates')
    rows = []
    for cohort, candidates in wall.items():
        for combo, count, pct in candidates:
            rows.append((cohort, json.dumps(combo), count, pct))
    c.executemany('INSERT INTO wall_candidates VALUES (?,?,?,?)', rows)

    conn.commit()
    conn.close()
    print(f"Guardado: {db_path}", flush=True)
    log.info(f"SQLite guardado: {db_path}")


if __name__ == '__main__':
    print("=== PIPELINE SNT — TCGA CORPUS ===\n", flush=True)
    patient_results, baselines = run_pipeline()
    cooc = build_cooccurrence_matrix(patient_results)
    wall = find_five_event_wall(patient_results)
    summary = generate_summary(patient_results, cooc, wall)
    save_to_sqlite(patient_results, wall)

    print("\n=== RESUMEN FINAL ===")
    for cohort, stats in summary.items():
        print(f"[{cohort}] n={stats['n_patients']} | "
              f"anom={stats['pct_with_anomaly']}% | "
              f"5+eventos={stats['pct_multi_event']}% | "
              f"wall_candidates={stats['n_wall_candidates']}")

    print("\n✓ Pipeline completado.")
