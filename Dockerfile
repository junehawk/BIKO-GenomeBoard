FROM python:3.12-slim

# WeasyPrint system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgdk-pixbuf2.0-0 \
    libffi-dev \
    libcairo2 \
    fonts-noto-cjk \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY config.yaml .
COPY scripts/ scripts/
COPY data/ data/
COPY templates/ templates/

# Create output and cache directories
RUN mkdir -p output data/db data/cache

# Volume mounts for:
# - /app/data/db: local ClinVar/gnomAD SQLite databases
# - /app/data/cache: variant response cache
# - /app/output: generated reports
# - /app/input: input VCF files
VOLUME ["/app/data/db", "/app/data/cache", "/app/output", "/app/input"]

# Default: show help
ENTRYPOINT ["python", "scripts/orchestrate.py"]
CMD ["--help"]
