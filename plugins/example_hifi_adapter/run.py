#!/usr/bin/env python3
"""Development-only HiFi adapter proving the long-read external contract."""

from __future__ import annotations

import sys
from pathlib import Path


# This adapter deliberately reuses only the deterministic fixture writer.  It
# is not a scientific long-read caller and is never listed as a formal tool.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "example_genotyper"))
from run import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
