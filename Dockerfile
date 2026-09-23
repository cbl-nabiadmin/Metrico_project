# Base Image: Lightweight Debian Linux pre-installed with Python 3.10
FROM python:3.10-slim

# Metadata
LABEL maintainer="Ardhendu"
LABEL description="Metrico: Containerized Functional Metagenomics Visualization and Statistical Analytics Pipeline"

# Prevent Python from writing .pyc files and enable unbuffered logging
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install system-level compilers required for C-extensions (scikit-bio, scipy)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgomp1 \
    && rm -rf /var/lib/apt-get/lists/*

# Set working directory inside the container
WORKDIR /Metrico_project

# Copy package requirements first to leverage Docker layer caching
COPY requirements.txt .

# Install Python packages without storing local caches
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy all scripts, licenses, default input files, and directories into container
COPY . /Metrico_project

# Set execution permissions on the master runner
RUN chmod +x /Metrico_project/run_all.py

# Default entry execution command
CMD ["python3", "/Metrico_project/run_all.py"]
