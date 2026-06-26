"""
mock_main.py — SNT Mock Integration Services
============================================
Simulates hospital SRE integrations:
  POST /jira/create_ticket
  POST /slack/notify_team
  POST /email/notify_reporter
  GET  /health

Each endpoint logs the incoming payload, simulates network latency,
and returns a synthetic response ID.

Author  : SNT Genomic Analyzer Team
License : MIT
"""

from __future__ import annotations

import logging
import sys
import time
import random
import string
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import uvicorn
from fastapi import FastAPI, Request, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# ── Logging ──────────────────────────────────────────────────────────────────
import os as _os
_log_dir  = _os.path.dirname(_os.path.abspath(__file__))
_log_file = "/data/mock_services.log" if _os.path.exists("/data") and _os.access("/data", _os.W_OK) \
            else _os.path.join(_log_dir, "mock_services.log")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(_log_file, mode="a"),
    ],
)
logger = logging.getLogger("SNT.MockServices")

# ── FastAPI app ───────────────────────────────────────────────────────────────
app = FastAPI(
    title="SNT Mock Integration Services",
    description=(
        "Simulated hospital/SRE integrations for the SNT Genomic Analyzer. "
        "Returns synthetic IDs and logs every call for observability."
    ),
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Request/Response Models ───────────────────────────────────────────────────

class JiraTicketRequest(BaseModel):
    project:     str
    summary:     str
    description: Optional[str] = ""
    priority:    Optional[str] = "Medium"
    labels:      Optional[list[str]] = []


class JiraTicketResponse(BaseModel):
    id:         str
    key:        str
    status:     str
    url:        str
    created_at: str


class SlackNotifyRequest(BaseModel):
    channel: str
    text:    str
    blocks:  Optional[list[dict[str, Any]]] = None


class SlackNotifyResponse(BaseModel):
    message_id: str
    channel:    str
    ok:         bool
    ts:         str


class EmailNotifyRequest(BaseModel):
    to:      str
    subject: str
    body:    str
    cc:      Optional[list[str]] = []


class EmailNotifyResponse(BaseModel):
    id:          str
    to:          str
    subject:     str
    queued_at:   str
    status:      str


# ── Helpers ───────────────────────────────────────────────────────────────────

def _random_id(prefix: str, length: int = 6) -> str:
    suffix = "".join(random.choices(string.digits, k=length))
    return f"{prefix}-{suffix}"


def _log_request(service: str, payload: Any) -> None:
    logger.info("[%s] Incoming request: %s", service, payload)


def _simulate_latency(service: str, min_ms: int = 80, max_ms: int = 400) -> None:
    delay = random.randint(min_ms, max_ms) / 1000
    logger.debug("[%s] Simulating network latency: %.0f ms", service, delay * 1000)
    time.sleep(delay)


# ── Middleware — Request logging ──────────────────────────────────────────────

@app.middleware("http")
async def log_all_requests(request: Request, call_next):
    logger.debug(
        "[HTTP] %s %s | Headers: %s",
        request.method,
        request.url.path,
        dict(request.headers),
    )
    response = await call_next(request)
    logger.debug("[HTTP] Response status: %d", response.status_code)
    return response


# ── Endpoints ─────────────────────────────────────────────────────────────────

@app.get("/health")
def health_check() -> dict[str, str]:
    """Health probe — used by Docker healthcheck and the agent."""
    logger.debug("[HEALTH] Health check ping received.")
    return {"status": "ok", "service": "snt-mock-services", "time": datetime.utcnow().isoformat()}


@app.post("/jira/create_ticket", response_model=JiraTicketResponse)
def create_jira_ticket(payload: JiraTicketRequest) -> JiraTicketResponse:
    """
    Simulate Jira ticket creation.
    Logs the full payload and returns a synthetic issue key.
    """
    _log_request("JIRA", payload.dict())
    _simulate_latency("JIRA", 150, 500)

    issue_id  = _random_id("GENOMICS", 5)
    issue_key = f"GENOMICS-{random.randint(1000, 9999)}"

    logger.info(
        "[JIRA] ✅ Ticket created | Key=%s | Priority=%s | Summary=%s",
        issue_key, payload.priority, payload.summary[:80],
    )

    return JiraTicketResponse(
        id=issue_id,
        key=issue_key,
        status="Open",
        url=f"https://jira.hospital.local/browse/{issue_key}",
        created_at=datetime.utcnow().isoformat(),
    )


@app.post("/slack/notify_team", response_model=SlackNotifyResponse)
def notify_slack_team(payload: SlackNotifyRequest) -> SlackNotifyResponse:
    """
    Simulate Slack message delivery to an oncology channel.
    """
    _log_request("SLACK", payload.dict())
    _simulate_latency("SLACK", 80, 300)

    ts = f"{int(time.time())}.{random.randint(100000, 999999)}"
    msg_id = str(uuid.uuid4())[:8].upper()

    logger.info(
        "[SLACK] ✅ Message delivered | Channel=%s | MsgID=%s | Text=%s",
        payload.channel, msg_id, payload.text[:100],
    )

    return SlackNotifyResponse(
        message_id=msg_id,
        channel=payload.channel,
        ok=True,
        ts=ts,
    )


@app.post("/email/notify_reporter", response_model=EmailNotifyResponse)
def notify_reporter_email(payload: EmailNotifyRequest) -> EmailNotifyResponse:
    """
    Simulate email dispatch to the oncology reporting team.
    """
    _log_request("EMAIL", payload.dict())
    _simulate_latency("EMAIL", 100, 400)

    email_id = f"MSG-{uuid.uuid4().hex[:12].upper()}"

    logger.info(
        "[EMAIL] ✅ Email queued | To=%s | Subject=%s | ID=%s",
        payload.to, payload.subject, email_id,
    )

    return EmailNotifyResponse(
        id=email_id,
        to=payload.to,
        subject=payload.subject,
        queued_at=datetime.utcnow().isoformat(),
        status="queued",
    )


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    logger.info("=" * 60)
    logger.info("SNT Mock Services starting on port 8081...")
    logger.info("=" * 60)
    uvicorn.run(app, host="0.0.0.0", port=8081, log_level="debug")
