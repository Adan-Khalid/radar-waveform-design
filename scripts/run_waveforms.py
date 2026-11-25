"""Pulse compression metrics of every waveform at both design points.

Writes reports/waveforms/metrics.csv, snr_gain.csv and ofdm_spread.csv.
"""

from __future__ import annotations

import numpy as np

from radarwave import waveforms as wf
from radarwave.config import REPORTS, load_config, write_csv
from radarwave.metrics import pulse_metrics

OUT = REPORTS / "waveforms"
C = 299_792_458.0


def snr_gain_db(s: np.ndarray, rng, trials: int = 20000, amp: float = 0.3) -> float:
    """Measured output SNR over input SNR per sample, matched filter, unit noise power."""
    n = s.size
    w = (rng.standard_normal((trials, n)) + 1j * rng.standard_normal((trials, n))) / np.sqrt(2)
    y_noise = w @ np.conj(s)
    out_snr = amp**2 * np.abs(np.sum(np.abs(s) ** 2)) ** 2 / np.var(y_noise)
    in_snr = amp**2 * np.mean(np.abs(s) ** 2)
    return float(10 * np.log10(out_snr / in_snr))


def main():
    cfg = load_config("waveforms")
    spc = cfg["analysis"]["samples_per_cell"]
    seed = cfg["analysis"]["seed"]
    b_hz = cfg["example_radar"]["bandwidth_hz"]
    cell_m = C / (2 * b_hz)
    rng = np.random.default_rng(seed)

    rows, gains = [], []
    for point, dp in cfg["design_points"].items():
        bt = dp["time_bandwidth"]
        for spec in dp["waveforms"]:
            s = wf.make_waveform(spec, bt, spc)
            m = pulse_metrics(s, samples_per_cell=spc)
            s1 = wf.make_waveform(spec, bt, 1)
            rows.append(
                {
                    "design_point": point,
                    "time_bandwidth": bt,
                    "waveform": wf.label(spec, bt),
                    "samples_per_cell": spc,
                    "psl_db": m["psl_db"],
                    "isl_db": m["isl_db"],
                    "width_3db_cells": m["width_3db_cells"],
                    "resolution_3db_m": m["width_3db_cells"] * cell_m,
                    "papr_db": wf.papr_db(s),
                    "papr_sampled_at_b_db": wf.papr_db(s1),
                }
            )
            gains.append(
                {
                    "design_point": point,
                    "waveform": wf.label(spec, bt),
                    "samples": s1.size,
                    "snr_gain_db": snr_gain_db(s1, rng),
                    "theory_db": 10 * np.log10(bt),
                }
            )
    write_csv(OUT / "metrics.csv", rows)
    write_csv(OUT / "snr_gain.csv", gains)

    # OFDM sidelobes depend on the random symbols: spread over many draws.
    draws = cfg["analysis"]["ofdm_draws"]
    spread = []
    for point, dp in cfg["design_points"].items():
        bt = dp["time_bandwidth"]
        vals = []
        for k in range(draws):
            s = wf.ofdm(wf.qpsk_symbols(bt, seed + 1000 + k), spc)
            m = pulse_metrics(s, samples_per_cell=spc)
            s1 = wf.ofdm(wf.qpsk_symbols(bt, seed + 1000 + k), 1)
            vals.append((m["psl_db"], m["isl_db"], wf.papr_db(s), wf.papr_db(s1)))
        v = np.array(vals)
        for j, name in enumerate(["psl_db", "isl_db", "papr_db", "papr_sampled_at_b_db"]):
            p5, p50, p95 = np.percentile(v[:, j], [5, 50, 95])
            spread.append(
                {
                    "design_point": point,
                    "subcarriers": bt,
                    "metric": name,
                    "draws": draws,
                    "p5": p5,
                    "median": p50,
                    "p95": p95,
                }
            )
    write_csv(OUT / "ofdm_spread.csv", spread)

    for r in rows:
        print(
            f"{r['waveform']:>14}: PSL {r['psl_db']:7.2f} dB  ISL {r['isl_db']:7.2f} dB  "
            f"width {r['width_3db_cells']:.3f} cells  PAPR {r['papr_db']:.2f} dB"
        )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
