#!/usr/bin/env python3
"""Explicit entry point for the separately frozen official joint-noise policy."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run import main

if __name__ == "__main__":
    raise SystemExit(main(noise_genotyping=True))
