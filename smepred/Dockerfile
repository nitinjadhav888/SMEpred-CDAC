# HelixZero-CMS — Docker image for ICE Cloud / any container platform
# Uses Python 3.11 slim for small image size

FROM python:3.11-slim

WORKDIR /app

# Install system lib for LightGBM (OpenMP)
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# Copy dependency manifest first (leverages Docker layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code and data
# We copy only what's needed for serving the API
COPY api/ api/
COPY src/ src/
COPY cli/ cli/
COPY data/modification_codes.json data/
COPY data/oligoformer/ data/oligoformer/
COPY models/ models/
COPY app.html .

# Expose the FastAPI default port
EXPOSE 8000

# Start the server
# ICE Cloud / Railway will inject $PORT env variable
CMD uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8000}
