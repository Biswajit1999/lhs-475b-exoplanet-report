"""Verify the repository copies of the Zenodo products."""

import hashlib
import json
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"


def main() -> None:
    manifest = json.loads((DATA / "source_manifest.json").read_text(encoding="utf-8"))
    failures = []
    for relative, expected in manifest["files"].items():
        path = DATA / relative
        actual = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        if actual != expected:
            failures.append(f"{relative}: expected {expected}, got {actual}")
    if failures:
        raise SystemExit("Source validation failed:\n" + "\n".join(failures))
    print(f"Verified {len(manifest['files'])} files from {manifest['record']}")


if __name__ == "__main__":
    main()
