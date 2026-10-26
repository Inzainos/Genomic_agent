# CHANGELOG — SNT Genomic Topologic Analyzer

Cronológico, inverso. Todo lo listado fue **medido corriendo**, no leído.

---

## 2026-10-25 — Alineación con SNT v2.7.1 (tela de incertidumbre) + 3 bugs

### Alineación con el marco (SNT v2.7.1 / marco v34)
- **Fin de la "confianza" inventada** (regla 8): los `confidence_score`
  hardcodeados (0.92, 0.95…) ya NO se muestran como "Conf: 92%" al LLM ni
  en la UI. Eran valores decorativos sin calibración de cohortes — el
  mismo patrón que Titan retiró ("asertividad 99.9%").
- **`evidence_tier` nuevo** (migración forward-only): `EMPIRICAL_5EVENT`
  (19 firmas del muro de 5 eventos, n=2,746 pacientes TCGA) vs
  `CURATED_LITERATURE` (25 firmas curadas de literatura). El triage
  reporta Z medido + tier; el ranking de la UI ordena por |Z|.
- **Marco de medición θ en el prompt LLM** (Axioma 0.2): B (razón R por
  par), V (proxy declarado m=R; u=NULL por diseño), ρ (1 sección
  transversal), W (un punto temporal), λ (|Z|>2.5 contra baseline n=40).
  El prompt declara que el diagnóstico es una sombra Π_θ — distinta
  patología puede producir la misma sombra.
- **φ honesto**: el índice de fricción F se etiqueta en el prompt como
  "φ proxy, umbral de diseño 0.5, no calibrado" (marco v34: φ postulado,
  no medido).
- README: ranking por |Z| + tier de evidencia, sin porcentajes inventados.

### Bugs corregidos (todos medidos antes del fix)
1. **Firmas duplicadas 4×** en la DB del volumen (176 filas / 44
   distintas): `INSERT OR IGNORE` sin restricción UNIQUE nunca ignoraba —
   cada `compose up` re-sembraba. Fix: dedup (id mínimo) + índice único
   `ux_disease_signature`. Verificado idempotente: 44/44 tras rebuild.
   (timeseries auditada por el mismo defecto: 20/20 correcta — 5 pares ×
   4 timepoints, sin bug.)
2. **DNS embebido de Docker roto en WSL2** (apt exit 100 en ambos
   builds; medido: `getent hosts deb.debian.org` FAIL dentro del
   contenedor, OK con `--dns 8.8.8.8` y OK con `--network=host`): builds
   del compose con `network: host` documentado en el propio compose.
3. **`.env.fractal` sin gitignore** (archivo cifrado FRCTL3, 600): el
   patrón `*.env` no matchea. Añadido `.env.fractal`/`.env.*` al
   `.gitignore` (verificado con `git check-ignore`).

### Verificación punta a punta
- DB construida local + en volumen: migración tier + dedup ✓ (19
  EMPIRICAL / 25 CURATED / 44 totales, idempotente).
- Triage smoke: expresión sana (medias del baseline) → ~0-3
  confirmaciones de banda; perturbada (CDK4 12×) → TNBC confirmada con
  tier y Z.
- Prompt LLM verificado programáticamente: `Conf:` ELIMINADO, marco θ
  PRESENTE, `u = NULL` PRESENTE, φ proxy PRESENTE (con resultado ACO).
- Servicios: `snt_agent_core` y `snt_mock_services` healthy, loopback
  127.0.0.1:8501/:8081, HTTP 200 ambos.

---

## 2026-10-20 — Separación como proyecto independiente

- Migrado de `The-shadow-Node-Theory/genomic_agent/` con 45 commits
  preservados (`git subtree split`); docker compose solo loopback;
  alineado con release 2.5.2 / marco v34 / auditoría v32.
