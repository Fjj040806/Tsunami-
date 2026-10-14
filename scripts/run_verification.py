"""Portable numerical test entry point; no editable install required."""
from pathlib import Path
import sys
import unittest
import json
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

if __name__ == "__main__":
    suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    output = ROOT / "results"
    output.mkdir(exist_ok=True)
    (output / "test_summary.json").write_text(json.dumps({
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "tests_run": result.testsRun, "failures": len(result.failures),
        "errors": len(result.errors), "successful": result.wasSuccessful(),
        "failure_details": [text for _, text in result.failures + result.errors],
        "command": "python scripts/run_verification.py",
    }, indent=2), encoding="utf-8")
    sys.exit(0 if result.wasSuccessful() else 1)
