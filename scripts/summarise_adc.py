"""Derived tables from the ADC runs (no new simulation).

1. reports/adc/a2_strong_echo.csv: the A2 losses read as "a strong echo elsewhere in the
   receive window sets the ADC full scale". The gain is set so the strong echo's largest
   complex sample magnitude just reaches full scale (an upper bound on its I or Q peak):
   A / sigma_c = sqrt(2 * SNR_strong * PAPR), with SNR_strong the strong echo's mean SNR
   per sample. A constant-envelope echo has PAPR 1;
   the OFDM 64 pulse has the PAPR of its samples (reports/waveforms/metrics.csv). Losses are
   interpolated linearly in dB between the A2 grid points; outside the grid they are blank.
2. reports/adc/a1_replicate_check.csv: A1 (seed 2026) against the replicate (seed 2027)
   for LFM BT 64, receiver by receiver.
"""

from __future__ import annotations

import numpy as np

from radarwave.config import REPORTS, read_csv, write_csv

STRONG_SNR_DB = [0, 10, 20, 30, 40]


def _interp(x, y, v):
    """Linear interpolation between the two grid points around v; blank if either is
    missing (outside the grid, or a receiver that never reached Pd 0.9)."""
    if not x[0] <= v <= x[-1]:
        return ""
    j = min(int(np.searchsorted(x, v, side="right")) - 1, x.size - 2)
    if not (np.isfinite(y[j]) and np.isfinite(y[j + 1])):
        return ""
    return float(y[j] + (v - x[j]) * (y[j + 1] - y[j]) / (x[j + 1] - x[j]))


def strong_echo_table():
    rows = read_csv(REPORTS / "adc/a2_losses.csv")
    metrics = read_csv(REPORTS / "waveforms/metrics.csv")
    papr = {
        "constant envelope (LFM, phase codes)": 0.0,
        "OFDM 64": float(
            next(m["papr_sampled_at_b_db"] for m in metrics if m["waveform"] == "OFDM 64")
        ),
    }
    out = []
    groups = {}
    for r in rows:
        if r["bits"] == "":
            continue
        key = (int(r["bits"]), r["kind"], r["dither"])
        groups.setdefault(key, []).append((float(r["a_over_sigma_db"]), float(r["loss_db"])))
    for (bits, kind, dither), pts in sorted(groups.items()):
        x, y = np.array(sorted(pts)).T
        for echo, p in papr.items():
            row = {"bits": bits, "kind": kind, "dither": dither, "strong_echo": echo, "papr_db": p}
            for snr in STRONG_SNR_DB:
                a_db = snr + 10 * np.log10(2) + p
                row[f"loss_db_at_{snr}db"] = _interp(x, y, a_db)
            out.append(row)
    write_csv(REPORTS / "adc/a2_strong_echo.csv", out)
    print(f"wrote {REPORTS / 'adc/a2_strong_echo.csv'}")


def replicate_table():
    a = {
        r["receiver"]: r
        for r in read_csv(REPORTS / "adc/a1_losses.csv")
        if r["waveform"] == "LFM (BT 64)"
    }
    b = {r["receiver"]: r for r in read_csv(REPORTS / "adc/a1rep_losses.csv")}
    out = []
    for rx, ra in a.items():
        rb = b[rx]
        la, lb = float(ra["loss_db"]), float(rb["loss_db"])
        half = np.hypot(
            (float(ra["loss_hi"]) - float(ra["loss_lo"])) / 2,
            (float(rb["loss_hi"]) - float(rb["loss_lo"])) / 2,
        )
        out.append(
            {
                "receiver": rx,
                "loss_seed_2026": la,
                "loss_seed_2027": lb,
                "difference": lb - la,
                "combined_half_width": half,
                "inside": bool(abs(lb - la) <= half),
                "required_snr_seed_2026": float(ra["required_snr_db"]),
                "required_snr_seed_2027": float(rb["required_snr_db"]),
            }
        )
    write_csv(REPORTS / "adc/a1_replicate_check.csv", out)
    inside = sum(r["inside"] for r in out)
    print(f"replicate: {inside} of {len(out)} differences inside the combined intervals")


def main():
    if (REPORTS / "adc/a2_losses.csv").exists():
        strong_echo_table()
    if (REPORTS / "adc/a1rep_losses.csv").exists():
        replicate_table()


if __name__ == "__main__":
    main()
