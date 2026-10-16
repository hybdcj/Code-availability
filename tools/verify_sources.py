"""Verify byte-for-byte source provenance; does not execute analyses."""
from pathlib import Path
import csv
import hashlib

root = Path(__file__).resolve().parents[1]
with (root / "source-manifest.csv").open(encoding="utf-8", newline="") as fh:
    rows = list(csv.DictReader(fh))
errors = []
for row in rows:
    path = (root / row["repository_path"]).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        errors.append(row["repository_path"] + ": missing/invalid path")
        continue
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != row["sha256"] or len(payload) != int(row["bytes"]):
        errors.append(row["repository_path"] + ": source differs from manifest")
if errors:
    raise SystemExit("\n".join(errors))
print(f"Verified {len(rows)} source files. No analyses executed.")
