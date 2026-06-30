# AGENTS_USE.md — SNT Agent: Use Cases, Security & Observability

## What This Agent Does

The SNT Genomic Topologic Analyzer Agent is an **autonomous multi-step reasoning pipeline** that converts raw RNA-seq expression data into a clinical genomic report without human intervention in the analysis loop.

The agent does not rely on a chat-loop. It executes a deterministic pipeline:

```
Input (text + optional CSV)
  → Guardrails Verification
  → Expression Data Loading / Ingestion
  → Level 1 Clinical Triage         (fast O(K) oracle cross-reference)
  → Level 2 Deep Block Scanner      (chromosome-by-chromosome Z-sweep)
  → LLM Diagnosis Generation        (Claude 3.5 Sonnet via OpenRouter)
  → Notification Dispatch           (Jira + Slack + Email + Telegram)
  → Structured Report Output
```

---

## Agent Use Cases

### Use Case 1 — Hospital Oncology Triage

**Scenario**: A patient presents with a suspicious breast mass. The pathology lab has produced an RNA-seq panel. The oncologist uploads the CSV and describes the clinical presentation.

**Agent behaviour**:
1. Level 1 cross-references against all 9 known disease signatures in milliseconds
2. Detects `MYC→CDK4 Z=+4.82` (LEAPFROG), `TP53→CDKN1A Z=-3.91` (HUB_COLLAPSE) — confirming TNBC signature
3. Dispatches Jira ticket to the genomics operations board
4. Sends Slack alert to `#oncology-alerts`
5. Delivers full LLM-generated medical report to the oncologist's email
6. Pushes a real Telegram push notification to the senior oncologist's mobile

**Value**: Turnaround from upload to report: under 60 seconds vs. days for manual genomic interpretation.

### Use Case 2 — Research Lab Biomarker Discovery

**Scenario**: A research lab has sequenced 50 patient tumors with no known diagnosis. They want to discover novel transcriptomic patterns.

**Agent behaviour**:
1. Level 1 returns no strong confirmed matches (or partial matches)
2. Level 2 scans all 10 chromosome blocks
3. Detects 6 orphan anomalies — hub-satellite decouplings with Z > 2.5 that do not match any known disease signature
4. LLM interprets these as potential novel biomarkers with suggested follow-up experiments

**Value**: The agent discovers what existing databases cannot tell you.

### Use Case 3 — Longitudinal Patient Monitoring

**Scenario**: The same patient's RNA-seq is analysed at 3-month intervals to track treatment response.

**Agent behaviour**:
- Each analysis run produces a new report with a new patient_id (e.g., `PX-001-T0`, `PX-001-T3`, `PX-001-T6`)
- The Z-Score trajectories across runs indicate whether anomalous hubs are normalising (treatment working) or worsening
- Each run dispatches an automated Jira ticket with the delta summary

---

## Security Architecture — Guardrails

### Threat Model

The agent accepts free-text clinical notes from hospital staff, patients, or research teams. This creates an attack surface for **prompt injection** — attempts to rewrite the agent's medical instructions by embedding adversarial commands in the input.

### Defence Layer: `verify_guardrails()`

Implemented in `agent_logic.py` before any LLM call:

**Stage 1 — Structural Validation**
- Empty input rejected
- Inputs over 4,000 characters rejected (prevents context overflow attacks)

**Stage 2 — Pattern Blacklist (16 signatures)**

| Pattern | Attack Type |
|---|---|
| `ignore previous instructions` | Classic instruction override |
| `ignore all instructions` | Variation override |
| `jailbreak` | Capability bypass |
| `act as` | Role impersonation |
| `you are now` | Identity hijacking |
| `forget your instructions` | Memory wipe attempt |
| `new instructions:` | Instruction injection |
| `system prompt:` | System context extraction |
| `bypass` / `override` | Access control circumvention |
| `reveal your prompt` | Prompt extraction |
| `print your system` | System context leakage |
| `sudo` | Privilege escalation |
| `[[` / `]]` | Template injection delimiters |
| `<\|im_start\|>` / `<\|system\|>` | Token-level injection |

**Stage 3 — Structural Wrapping**

Even if a malicious string passes Stage 2, the user input is never passed directly as a system prompt. It is embedded inside a structured medical prompt template:

```python
f"""You are a senior clinical geneticist...
PATIENT: {patient_id}
CLINICAL NOTES: {clinical_notes}   # ← user input, sandboxed here
=== SNT TWO-LEVEL ANALYSIS RESULTS ===
{analysis_data}                    # ← only our structured data after this
"""
```

The LLM sees the user input as a *data field*, not as instructions.

---

## Observability — Log Pipeline

All services emit structured logs to both stdout and `/data/*.log`.

### Log Files

| Service | Log File | Rotation |
|---|---|---|
| `db_builder` | `/data/db_builder.log` | Single run |
| `mock_services` | `/data/mock_services.log` | Append |
| `agent_core` | `/data/agent_core.log` | Append |
| Streamlit UI | `/data/streamlit_ui.log` | Append |

### Log Event Taxonomy

Every key pipeline transition has a tagged log event:

```
[GUARDRAILS]  Input inspection result
[DB]          Database read/write operations
[LEVEL-1]     Triage cross-reference per signature
[LEVEL-2]     Block scanner per chromosome
[LLM]         OpenRouter request/response lifecycle
[NOTIFY]      Mock service call results
[TELEGRAM]    Real notification delivery
[HTTP]        All HTTP requests through mock services
[PIPELINE]    Top-level stage transitions with timing
```

### Example Log Trace (single analysis run)

```
10:23:38 | INFO  | [GUARDRAILS] Input cleared. Proceeding to analysis.
10:23:38 | INFO  | [DB] Loaded 54 gene expression values.
10:23:38 | INFO  | [PIPELINE] ── LEVEL 1 TRIAGE ──────────────────────────
10:23:38 | INFO  | [LEVEL-1] chr8  | Breast_Cancer_Basal_TNBC | Hub=MYC   Sat=CDK4    Z=+4.821 ✓ CONFIRMED
10:23:38 | INFO  | [LEVEL-1] chr17 | Breast_Cancer_Basal_TNBC | Hub=TP53  Sat=CDKN1A  Z=-3.914 ✓ CONFIRMED
10:23:38 | INFO  | [LEVEL-1] Triage complete. 8/25 signatures confirmed.
10:23:38 | INFO  | [PIPELINE] ── LEVEL 2 DEEP SCANNER ────────────────────
10:23:38 | INFO  | [LEVEL-2] ▶ Scanning block: chr1
10:23:38 | INFO  | [LEVEL-2] ✓ Block chr1 complete. 2 pairs scanned, 0 anomalies.
10:23:38 | INFO  | [LEVEL-2] ▶ Scanning block: chr3
10:23:38 | WARNING| [LEVEL-2] ORPHAN ANOMALY | chr3 | VHL→HIF1A | Z=+2.87 | Type=LEAPFROG
10:23:39 | INFO  | [PIPELINE] ── LLM DIAGNOSIS ──────────────────────────
10:23:42 | INFO  | [LLM] Diagnosis received (1842 chars).
10:23:42 | INFO  | [JIRA] ✅ Ticket created | Key=GENOMICS-4821
10:23:42 | INFO  | [SLACK] ✅ Message delivered | Channel=#oncology-alerts
10:23:42 | INFO  | [EMAIL] ✅ Email queued | To=oncology-team@hospital.local
10:23:42 | INFO  | [TELEGRAM] Message sent successfully.
10:23:42 | INFO  | [PIPELINE] Analysis complete for DEMO-PX-001 in 4.31 seconds.
```

### Viewing Logs

```bash
# All services
docker compose logs -f

# Specific service
docker compose logs -f agent_core
docker compose logs -f mock_services

# Inside the UI
# Navigate to the "System Logs" tab in Streamlit
```

---

## Multimodal Input Specification

The Streamlit UI accepts two input modalities simultaneously:

### Modality 1 — Text (Clinical Notes)
- Free-text description of clinical presentation, symptoms, family history
- Processed through Guardrails before any LLM interaction
- Max 4,000 characters

### Modality 2 — CSV File (RNA-seq Expression)
- Format: two columns, `gene_id` and `tpm_value`
- Header row is optional (auto-detected)
- Gene IDs are normalised to uppercase on ingest
- Malformed rows are skipped with a WARNING log and reported to the user
- If no file is uploaded, the system falls back to the patient's pre-loaded data in `patient_expression`

```csv
gene_id,tpm_value
MYC,450.3
TP53,22.1
BRCA1,30.4
CDK4,280.1
...
```


---

## ACO-A & 5-Event Wall (SNT v30)

El agente integra el framework de **Colapso Orbital Acoplado (ACO-A)**:

### Ejes ortogonales (b ⊥ Δ)
- **Eje b** (satelización): `R(t) = a·t^b` — dominancia hub vs satélite mientras corre el acoplamiento
- **Eje Δ** (absorción): `A(τ) = c·τ^Δ` — velocidad de absorción post-colapso del hub

### Fricción biológica F_bio
```
F_bio = mean(TPM_guardian_genes) / 100
Guardian genes: TP53, BRCA1, BRCA2, MLH1, ATM, CHEK2, RAD51, FANCD2, RB1, PTEN
F > 0.5  → alta fricción → Regulated_Orbital_Decay
F ≤ 0.5  → baja fricción → Catastrophic_Cliff / Cracquelure_Decay / Floor_Arrested
```

### Modos de colapso
| Modo | Condición | Descripción |
|------|-----------|-------------|
| Regulated_Orbital_Decay | F > 0.5 | Absorción suave, guardada por fricción |
| Cracquelure_Decay | F ≈ 0, gradual | Fragmentación progresiva |
| Floor_Arrested | F ≈ 0, abrupto, piso | Colapso abrupto con piso estable |
| Catastrophic_Cliff | F ≈ 0, abrupto, sin piso | Caída libre sin recuperación |
| Logistic_Sweep | magnitud acotada | Barrido logístico (ej. variante viral) |

### 5-Event Wall
Firma empírica de ≥5 colaps os hub simultáneos. Derivada de corpus TCGA n=2,746.
- Conecta con la Ley de Inevitabilidad: cuando ≥5 hubs colapsan, h(τ) → máximo
- Ver tab "🧬 5-Event Wall" en la UI Streamlit
- Endpoint REST: `GET /tcga/wall_summary`

### Uso CLI
```bash
# Modo clínico individual (ratio-based)
python run_ingesta.py --patient PX-001 --csv paciente.csv --mode ratio

# Modo poblacional (gene-based, TCGA-compatible)
python run_ingesta.py --patient PX-001 --csv paciente.csv --mode gene

# Export JSON con resultados ACO-A incluidos
python run_ingesta.py --patient PX-001 --notes "LUAD suspected" --json-out report.json
```
