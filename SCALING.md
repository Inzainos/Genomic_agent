# SCALING.md — SNT Genomic Analyzer: From 10 Patients to 10,000

## Current Architecture Limits

The current stack is designed for hackathon demonstration and single-institution clinical use. At this scale:

- SQLite handles concurrent reads well, but write contention begins above ~50 concurrent patients
- The Streamlit server handles one analysis thread at a time per worker process
- Mock services have no persistence — all notification history is in-memory logs
- The LLM call is synchronous and blocks the request thread

This document describes the evolution path from the current Docker Compose deployment to a production-grade hospital informatics platform.

---

## Phase 1 — Multi-Patient Concurrency (10 → 100 patients/day)

**Problem**: SQLite serialises writes. Multiple simultaneous uploads corrupt `patient_expression`.

**Solution**: Replace SQLite with PostgreSQL.

```yaml
# docker-compose.production.yml addition
  postgres:
    image: postgres:16
    environment:
      POSTGRES_DB: snt_genomic
      POSTGRES_USER: snt_agent
      POSTGRES_PASSWORD: ${POSTGRES_PASSWORD}
    volumes:
      - pg_data:/var/lib/postgresql/data
```

The `db_builder.py` schema is already parameterised — switching the connection string in `agent_logic.py` from `sqlite3` to `psycopg2` requires changing two lines.

**Problem**: Streamlit blocks on the LLM call (up to 30s per patient).

**Solution**: Introduce a task queue. Replace the synchronous `run_full_analysis()` call with a Celery task dispatched to a Redis broker.

```
User submits → FastAPI endpoint → Celery.delay(run_full_analysis) → Returns task_id
                                                   ↓
                                          Worker picks up task
                                                   ↓
                                     SSE stream result back to UI
```

---

## Phase 2 — High-Throughput Research (100 → 10,000 patients/day)

**Problem**: Level 2 Block Scanner is sequential per patient. At 10,000 patients × 10 chromosomes = 100,000 scan jobs/day.

**Solution**: Replace the for-loop with a distributed chromosome job queue.

```
Kafka Topic: snt.chromosome_scan_jobs
  → Partition by chromosome (chr1...chr22, chrX, chrY)
  → Consumer group: snt-level2-workers (autoscale 1–20 pods)
  → Results sink to: snt.orphan_anomalies (Kafka) → PostgreSQL
```

Each Kafka partition maps to one chromosome. The 22 autosomes + sex chromosomes = 24 natural partitions. Each worker only ever loads one chromosome's baseline data at a time — the RAM-safety property of the block scanner is preserved.

**Problem**: 10,000 × LLM calls/day = significant cost and latency.

**Solution**: Two-tier LLM strategy.

| Tier | Condition | Model |
|---|---|---|
| Fast Triage | Level 1 returns ≥1 confirmed match | `claude-haiku-3-5` (fast, cheap) |
| Deep Diagnosis | Novel orphan anomalies OR complex presentation | `claude-3.5-sonnet` (full capability) |

The routing logic lives in `agent_logic.py::call_llm_for_diagnosis()`. A simple conditional on `len(orphan_anomalies) > 0 or len(confirmed) == 0` selects the model.

---

## Phase 3 — National / Multi-Institutional Scale (10,000+ patients/day)

### Data Layer

Replace single PostgreSQL with a federated architecture:

```
Hospital A (EU)                 Hospital B (US)
  ├── Local patient_expression    ├── Local patient_expression
  └── Read replica baseline       └── Read replica baseline
                 ↓                              ↓
         Central baseline_network_reference (read-only)
         Central disease_snt_signatures      (read-only)
         Federated orphan_anomalies          (write from all)
```

The baseline and disease oracle tables are read-only reference data and can be cached at each institution. Only orphan anomalies discovered at each site need to propagate back to the central registry — this is the **scientific value accumulation loop**.

### Compute Layer

```
Kubernetes (GKE / EKS / AKS)
  ├── snt-db-service          PostgreSQL + PgBouncer (connection pooling)
  ├── snt-api (FastAPI)        Horizontal pod autoscaler: 2–50 replicas
  ├── snt-worker (Celery)      Autoscaler tied to Redis queue depth
  ├── snt-ui (Streamlit)       2–10 replicas behind nginx ingress
  └── snt-mock-services        Replace with real Jira/Slack webhooks
```

### Observability at Scale

Replace file-based logging with:

- **OpenTelemetry** spans for every pipeline stage
- **Prometheus** metrics: `snt_analyses_total`, `snt_anomalies_detected`, `snt_llm_latency_seconds`
- **Grafana** dashboard: real-time anomaly discovery rate, LLM cost tracker, per-chromosome scan throughput

The log format used in the current stack (`%(asctime)s | %(levelname)-8s | %(name)s | %(message)s`) is already structured and trivially convertible to JSON for ingestion into Loki or Elasticsearch.

---

## SNT-Specific Scaling Property

Unlike conventional diagnostic pipelines that scale linearly with genome size, the SNT Two-Level Architecture has a key **sub-linear scaling property**:

- Level 1 complexity is **O(K)** — it does not grow with genome size, only with the number of known diseases (currently 25 signatures across 9 diseases, a set that grows slowly)
- Level 2 complexity is **O(chromosomes × avg_pairs_per_chromosome)** — it is bounded and predictable; each chromosome block is processed independently and discarded

This means adding 100 new patients adds 100 Level 1 queries (fast) and 100 × 24 chromosome jobs (parallelisable). The system never loads the full genome for any single patient into RAM.

---

## Scaling Milestones Summary

| Scale | Bottleneck | Solution |
|---|---|---|
| 10 patients/day | None | Current Docker Compose |
| 100 patients/day | SQLite writes, Streamlit blocking | PostgreSQL + Celery/Redis |
| 1,000 patients/day | Sequential Level 2 scanner | Kafka-partitioned chromosome workers |
| 10,000 patients/day | LLM cost and latency | Model routing (Haiku vs Sonnet) |
| 100,000 patients/day | Single DB write path | Federated PostgreSQL + K8s HPA |
