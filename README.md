# 🧬 SNT Genomic Topologic Analyzer Agent

> **Functional genomic diagnosis through hub-satellite network topology — not sequence.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Docker](https://img.shields.io/badge/Docker-Compose-2496ED?logo=docker)](docker-compose.yml)
[![Python](https://img.shields.io/badge/Python-3.10-3776AB?logo=python)](https://python.org)
[![LLM](https://img.shields.io/badge/LLM-Claude%203.5%20Sonnet-6B48FF)](https://openrouter.ai)

---

## The Paradigm Shift

Traditional genomic diagnostics detect **structural mutations** — a wrong nucleotide, a deleted exon. The SNT Analyzer instead asks: *who stopped controlling whom?*

When a master regulator (hub) loses control of its downstream genes (satellites), that loss of control **is** the disease signal — often before a mutation is even detectable. This is **Functional Topological Diagnosis**.

| Dimension | Traditional Genomics | SNT Topology |
|---|---|---|
| Unit of analysis | Individual mutations (SNPs, CNVs) | Hub-satellite regulatory relationships |
| Signal type | Structural ("wrong letter") | Behavioural ("wrong control pattern") |
| Complexity | O(genome) per variant | O(K) for triage, O(chromosome) for scan |
| Novel discovery | Constrained to known variants | Orphan anomaly detection for new biomarkers |
| RAM footprint | Scales with genome width | Constant per chromosome block |

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    snt_network (Docker Bridge)                  │
│                                                                 │
│  ┌──────────────┐    ┌────────────────────┐    ┌────────────┐  │
│  │  db_builder  │───▶│    agent_core      │───▶│  /data     │  │
│  │  (SQLite)    │    │  Streamlit :8501   │    │  (volume)  │  │
│  └──────────────┘    │  agent_logic.py    │    └────────────┘  │
│                      │  app.py            │                     │
│  ┌──────────────┐    └────────┬───────────┘                     │
│  │mock_services │◀───────────┘                                  │
│  │  FastAPI     │    Jira / Slack / Email mocks                 │
│  │  :8081       │                                               │
│  └──────────────┘                                               │
└─────────────────────────────────────────────────────────────────┘
                              │
                    Real Telegram API ──▶ 📱 Oncologist
```

### Three Modules

**`/genomic_database`** — The SNT Storage Layer
- SQLite with 3 tables: healthy tissue reference network, patient RNA-seq expression, clinical disease oracle
- Seeded with biologically calibrated data (MYC regulon, TP53 regulon, BRCA1, EGFR, PIK3CA, and 8 more hubs across 10 chromosomes)
- 9 disease signatures covering TNBC, Lung Adenocarcinoma, Colorectal Cancer, Melanoma, GBM, Renal Cell Carcinoma, PDAC, Li-Fraumeni Syndrome, HBOC

**`/mock_services`** — Hospital Integration Layer (FastAPI :8081)
- `/jira/create_ticket` — Simulates Jira issue creation with realistic latency
- `/slack/notify_team` — Simulates Slack channel notification
- `/email/notify_reporter` — Simulates email dispatch to oncology team

**`/agent_core`** — Intelligence + UI (Streamlit :8501)
- `agent_logic.py` — Two-Level SNT Math Engine + LLM Orchestration
- `app.py` — Multimodal Streamlit interface (text + CSV upload)

---

## Two-Level Analysis Engine

### Level 1 — Clinical Triage (milliseconds)

Complexity: **O(K)** where K = number of known disease-gene pairs.

The engine loads only the disease oracle signatures and cross-references them against patient expression. Each hub-satellite pair is scored:

```
Z = (R_patient - μ_healthy) / σ_healthy
where R = TPM(satellite) / TPM(hub)
```

Pairs with |Z| > 2.5 are confirmed. Output: ranked list of matching diseases with confidence scores.

### Level 2 — Deep Block Scanner (background)

Iterates chromosome by chromosome. For each block:
1. Load hub-satellite pairs for that chromosome from the baseline
2. Compute Z-Scores for all pairs present in patient expression
3. Flag pairs with |Z| > 2.5 as **Orphan Anomalies** (not in any known disease)
4. Release block from RAM
5. Move to next chromosome

Output: novel topological anomalies — potential new biomarkers or undescribed disease subtypes.

---

## Three Anomaly Types

| Type | Biological Meaning | Z-Score Direction |
|---|---|---|
| **LEAPFROG** | Satellite decouples upward — activates independent pathways | Z >> +2.5 |
| **SATELLITE_CAPTURE** | Free gene drawn into hub dependency abnormally | Z << -2.5 |
| **HUB_COLLAPSE** | Master regulator loses control of its entire regulon | Hub TPM collapse |

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/Inzainos/The-shadow-Node-Theory.git
cd The-shadow-Node-Theory/genomic_agent

# 2. Configure
cp .env.example .env
# Edit .env with your OPENROUTER_API_KEY

# 3. Launch — everything runs in containers, zero host dependencies
docker compose --env-file .env up --build

# 4. Open
# Streamlit UI  →  http://localhost:8501
# Mock Services →  http://localhost:8081/docs
```

For the complete step-by-step guide including demo walkthrough and CSV format, see [QUICKGUIDE.md](QUICKGUIDE.md).

---

## Security — Guardrails

Every user input is validated before reaching the LLM:

- **Injection pattern blacklist** — 16 known prompt injection signatures neutralised
- **Length cap** — inputs over 4,000 characters are rejected
- **Empty input guard** — blank submissions are blocked at the API layer
- **No raw passthrough** — user text is wrapped in a structured medical prompt, never inserted as system instructions

See [AGENTS_USE.md](AGENTS_USE.md) for the full security model.

---

## Observability

All services log to both stdout (Docker logs) and shared `/data/*.log` files:

```
2025-01-15 10:23:41 | INFO     | SNT.AgentLogic | [LEVEL-1] ✓ CONFIRMED | chr8 | Breast_Cancer_Basal_TNBC | Hub=MYC     Sat=CDK4     Z=+4.821
2025-01-15 10:23:41 | WARNING  | SNT.AgentLogic | [LEVEL-2] ORPHAN ANOMALY | chr12 | KRAS→ERK1 | Z=+3.104 | Type=LEAPFROG
2025-01-15 10:23:42 | INFO     | SNT.AgentLogic | [LLM] Diagnosis received (1842 chars).
2025-01-15 10:23:42 | INFO     | SNT.AgentLogic | [JIRA] ✅ Ticket created | Key=GENOMICS-4821
2025-01-15 10:23:42 | INFO     | SNT.AgentLogic | [TELEGRAM] Message sent successfully.
```

Live logs are also visible in the **System Logs** tab of the Streamlit UI.

---

## Repository Structure

```
snt-genomic-analyzer/
├── docker-compose.yml          # Master orchestrator
├── .env.example                # Environment variable template
├── LICENSE                     # MIT
├── README.md
├── AGENTS_USE.md               # Agent use cases & security deep-dive
├── SCALING.md                  # Production scaling strategy
├── QUICKGUIDE.md               # Judge-facing 5-minute demo guide
├── genomic_database/
│   ├── Dockerfile
│   ├── requirements.txt
│   └── db_builder.py           # Schema + synthetic data seeder
├── mock_services/
│   ├── Dockerfile
│   ├── requirements.txt
│   └── mock_main.py            # FastAPI mock integrations
└── agent_core/
    ├── Dockerfile
    ├── requirements.txt
    ├── agent_logic.py          # SNT math engine + LLM + notifications
    └── app.py                  # Streamlit multimodal UI
```
