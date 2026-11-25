"""Settings files, report paths and small file helpers shared by the scripts."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIGS = REPO_ROOT / "configs"
REPORTS = REPO_ROOT / "reports"
STYLE = REPO_ROOT / "assets" / "style.mplstyle"


def load_config(name: str | Path) -> dict:
    """Load a YAML settings file by name ('adc') or by path."""
    path = Path(name)
    if path.suffix not in (".yaml", ".yml"):
        path = CONFIGS / f"{name}.yaml"
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def write_csv(path: str | Path, rows: list[dict], fields: list[str] | None = None, digits=8):
    """Write dict rows to CSV; floats rounded to `digits` significant digits."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k), digits) for k in fields})


def read_csv(path: str | Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _fmt(v, digits):
    if isinstance(v, (float,)) or hasattr(v, "dtype"):
        v = float(v)
        if math.isnan(v):
            return ""
        if math.isinf(v):
            return "inf" if v > 0 else "-inf"
        if v == 0.0:
            v = 0.0  # no "-0" in the tables
        return f"{v:.{digits}g}"
    return v


def use_style():
    """Matplotlib with the repository style, no display needed."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use(str(STYLE))
    return plt
