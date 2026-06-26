FROM python:3.10-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY mock_main.py .

RUN mkdir -p /data

EXPOSE 8081

CMD ["python", "mock_main.py"]
