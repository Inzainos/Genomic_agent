FROM python:3.10-slim

WORKDIR /app

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    sqlite3 \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies (cached layer)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Source files
COPY agent_logic.py .
COPY app.py .
COPY run_ingesta.py .
COPY data_sanitizer.py .

# Non-root user for security
RUN useradd -m -u 1000 sntuser && chown -R sntuser:sntuser /app
RUN mkdir -p /data && chown sntuser:sntuser /data
USER sntuser

EXPOSE 8501

HEALTHCHECK --interval=15s --timeout=10s --start-period=30s --retries=5 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8501/_stcore/health')" || exit 1

CMD ["streamlit", "run", "app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true", \
     "--browser.gatherUsageStats=false"]
