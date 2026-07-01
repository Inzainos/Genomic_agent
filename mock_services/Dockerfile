FROM python:3.10-slim

WORKDIR /app

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    socat \
    bubblewrap \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY mock_main.py .

# Non-root user for security
RUN useradd -m -u 1000 sntuser && chown -R sntuser:sntuser /app
RUN mkdir -p /data && chown sntuser:sntuser /data
USER sntuser

EXPOSE 8081

HEALTHCHECK --interval=10s --timeout=5s --start-period=10s --retries=5 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8081/health')" || exit 1

CMD ["python", "mock_main.py"]
