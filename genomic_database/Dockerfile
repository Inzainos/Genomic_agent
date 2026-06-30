FROM python:3.10-slim

WORKDIR /app

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    sqlite3 \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Source
COPY db_builder.py .
COPY hpa_db_builder.py .

# Non-root user for security
RUN useradd -m -u 1000 sntuser && chown -R sntuser:sntuser /app
RUN mkdir -p /data && chown sntuser:sntuser /data
USER sntuser

# Build the database on container start
CMD ["python", "db_builder.py"]
