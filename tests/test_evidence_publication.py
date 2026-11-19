"""Keep research evidence and its public page in sync with versioned inputs."""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_evidence_publication_is_current():
    snapshot = subprocess.run(
        ["git", "cat-file", "-e", "HEAD:evidence/manifest.json"],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    ) if (ROOT / ".git").exists() else None
    option = "--check-committed" if snapshot and snapshot.returncode == 0 else "--check"
    for script in ("scripts/evidence_catalog.py", "scripts/render_evidence_page.py"):
        result = subprocess.run(
            [sys.executable, script, option], cwd=ROOT, capture_output=True, text=True, check=False
        )
        assert result.returncode == 0, result.stdout + result.stderr
