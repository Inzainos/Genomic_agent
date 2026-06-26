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

# Data volume mountpoint
RUN mkdir -p /data

# Build the database on container start
CMD ["python", "db_builder.py"]
