#!/usr/bin/env python3
"""
ICE Cloud Deployment Guide for HelixZero-CMS
=============================================

This script does NOT deploy automatically — it prints step-by-step
instructions and optionally builds the Docker image locally.

Usage:
    python deploy_icecloud.py          # Show instructions only
    python deploy_icecloud.py --build  # Also build the Docker image
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path


STEPS = """
╔══════════════════════════════════════════════════════════════╗
║   HelixZero-CMS — ICE Cloud Deployment Steps               ║
║   C-DAC Pune | ICE Cloud → icecloud.in                     ║
╚══════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PREREQUISITES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. ICE Cloud account (you@cdac.in) — login at:
   → https://icecloud.in/me/dashboard

2. Docker Desktop installed on your machine:
   → https://www.docker.com/products/docker-desktop/

3. Open PowerShell / CMD as Administrator and run:
   → docker --version   (should show 24+)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 1: BUILD THE DOCKER IMAGE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  cd D:\\Helixx\\smepred
  docker build -t helixzero-cms:latest .

  ⏱ Takes ~2-3 minutes first time (downloads Python + installs deps)
  ✅ Final output: "Successfully tagged helixzero-cms:latest"

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 2: TEST LOCALLY (optional but recommended)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  docker run -d -p 8000:8000 --name helixzero helixzero-cms:latest

  Open → http://localhost:8000
  Test → curl http://localhost:8000/modifications

  Stop → docker stop helixzero && docker rm helixzero

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 3: TAG & PUSH TO ICE CLOUD REGISTRY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # Login to ICE Cloud container registry
  docker login icecloud.in  --username YOUR_CDAC_EMAIL

  # Tag the image with your ICE Cloud project namespace
  docker tag helixzero-cms:latest icecloud.in/YOUR_USERNAME/helixzero-cms:latest

  # Push the image
  docker push icecloud.in/YOUR_USERNAME/helixzero-cms:latest

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 4: DEPLOY ON ICE CLOUD DASHBOARD
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. Go to → https://icecloud.in/me/dashboard

  2. Click "Containers" or "Workloads" → "Create Deployment"

  3. Fill in:

     ┌─────────────────────────────────────────────┐
     │ Image:        icecloud.in/YOUR_USERNAME/    │
     │               helixzero-cms:latest           │
     │                                              │
     │ Name:         helixzero-cms                  │
     │                                              │
     │ Port:         8000                           │
     │                                              │
     │ Memory:       2048 MB   (2 GB recommended)   │
     │                                              │
     │ CPU:          1-2 vCPU                       │
     │                                              │
     │ Env Variable: PORT = 8000                    │
     └─────────────────────────────────────────────┘

  4. Click "Deploy"

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 5: VERIFY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  ICE Cloud will give you a URL like:
  → https://helixzero-cms-xxxx.icecloud.in

  Test it:
  → curl https://helixzero-cms-xxxx.icecloud.in/modifications

  If you see the modification list JSON → SUCCESS 🎉

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 6: USE THE APP
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Open the URL in browser → the app.html frontend loads.
  Paste a gene sequence → click "Rank" → get scores.

  API endpoints available:
  POST /rank          — Rank unmodified siRNA candidates
  POST /single-mod    — Single-modification scan
  POST /multi-mod-scan — Beam search (multi-mod)
  GET  /modifications — List all supported modifications

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DONE. 🚀
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""


def build_image():
    """Run `docker build` in the project root."""
    project_root = Path(__file__).parent
    print("🐳 Building Docker image: helixzero-cms:latest ...")
    result = subprocess.run(
        ["docker", "build", "-t", "helixzero-cms:latest", "."],
        cwd=project_root,
        capture_output=False,
    )
    if result.returncode == 0:
        print("\n✅ Build successful! Image: helixzero-cms:latest")
    else:
        print("\n❌ Build failed. Check Docker is running and try again.")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="ICE Cloud deployment helper")
    parser.add_argument("--build", action="store_true", help="Build Docker image")
    args = parser.parse_args()

    if args.build:
        build_image()

    print(STEPS)


if __name__ == "__main__":
    main()
