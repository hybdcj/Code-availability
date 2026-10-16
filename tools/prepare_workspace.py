"""Copy sources to an explicit local workspace; never execute or copy research data."""
from pathlib import Path
import argparse
import csv
import hashlib
import json

def prepare(config_path):
    repo = Path(__file__).resolve().parents[1]
    cfg = json.loads(Path(config_path).read_text(encoding="utf-8-sig"))
    keys = ("project_root", "scratch_root", "external_staging_root", "r_library")
    for key in keys:
        if not isinstance(cfg.get(key), str) or not Path(cfg[key]).is_absolute():
            raise ValueError(f"{key} must be an absolute path on this machine")
        cfg[key] = Path(cfg[key]).resolve().as_posix()
        if any(c in cfg[key] for c in "\"'\n\r"):
            raise ValueError(f"{key} cannot contain quotes or newlines")
    project = Path(cfg["project_root"])
    if project == repo or project.is_relative_to(repo):
        raise ValueError("Use a work directory outside the repository")
    mappings = {
        "D:/BaiduSyncdisk/科研/徐中豪": cfg["project_root"],
        "C:/y03tmp": cfg["scratch_root"],
        "D:/codex_tmp": cfg["external_staging_root"],
        "D:/R/library": cfg["r_library"],
    }
    with (repo / "source-manifest.csv").open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    planned = []
    for row in rows:
        source = (repo / row["repository_path"]).resolve()
        if not source.is_relative_to(repo):
            raise ValueError("Source path escapes repository")
        raw = source.read_bytes()
        if hashlib.sha256(raw).hexdigest() != row["sha256"]:
            raise ValueError(f"Source checksum mismatch: {source}")
        target = (project / row["project_relative_path"]).resolve()
        if not target.is_relative_to(project):
            raise ValueError("Destination path escapes workspace")
        text = raw.decode("utf-8-sig")
        for old, new in mappings.items():
            text = text.replace(old, new)
            text = text.replace(old.replace("/", "\\"), new)
        payload = text.encode("utf-8")
        if target.exists() and target.read_bytes() != payload:
            raise FileExistsError(f"Refusing to overwrite a different existing file: {target}")
        planned.append((target, payload, row["repository_path"]))
    # Validate the whole plan before writing any source file.
    for target, payload, _ in planned:
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            target.write_bytes(payload)
    for base, rels in ((project, ("output", "16.论文修订第三轮/_核查/输出")),
                       (Path(cfg["scratch_root"]), ("tmp", "out", "F34_delivery/TIFF", "F36_delivery/TIFF")),
                       (Path(cfg["external_staging_root"]), ("hrs", "elsa", "share", "ckm_stage"))):
        for rel in rels:
            (base / rel).mkdir(parents=True, exist_ok=True)
    report = [{"source": src, "target": path.as_posix(), "prepared_sha256": hashlib.sha256(payload).hexdigest()}
              for path, payload, src in planned]
    (project / "code-preparation-manifest.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Prepared {len(planned)} text files in {project}; no data copied or analysis executed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    prepare(args.config)
