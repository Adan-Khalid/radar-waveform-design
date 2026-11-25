"""Doppler tolerance of every waveform: peak loss, delay shift and PSL against Doppler.

Writes reports/ambiguity/doppler_tolerance.csv and doppler_cuts.csv.
"""

from __future__ import annotations

import numpy as np

from radarwave import waveforms as wf
from radarwave.ambiguity import doppler_cut_summary, doppler_tolerance, first_crossing
from radarwave.config import REPORTS, load_config, write_csv

OUT = REPORTS / "ambiguity"
C = 299_792_458.0
FINE_MAX = 4.0  # nu*T range of the fine scan


def main():
    cfg = load_config("waveforms")
    spc = cfg["analysis"]["samples_per_cell"]
    fc = cfg["example_radar"]["carrier_hz"]
    b_hz = cfg["example_radar"]["bandwidth_hz"]
    lam = C / fc

    summary, cuts = [], []
    for point, dp in cfg["design_points"].items():
        bt = dp["time_bandwidth"]
        t_pulse = bt / b_hz
        for spec in dp["waveforms"]:
            s = wf.make_waveform(spec, bt, spc)
            name = wf.label(spec, bt)
            summ, cut = doppler_tolerance(s, nu_t_max=FINE_MAX, step=0.01, samples_per_cell=spc)
            for key in ("nu_t_1db", "nu_t_3db"):
                if not np.isfinite(summ[key]):  # Doppler-tolerant: scan up to nu*T = BT
                    wide = np.round(np.arange(0, bt + 1e-9, 0.05), 10)
                    c2 = doppler_cut_summary(s, None, wide, spc)
                    level = -1.0 if key == "nu_t_1db" else -3.0
                    summ[key] = first_crossing(wide, c2["peak_loss_db"], level)
            row = {"design_point": point, "time_bandwidth": bt, "waveform": name, **summ}
            for key in ("nu_t_1db", "nu_t_3db"):
                nu = summ[key] / t_pulse
                row[key.replace("nu_t", "doppler_hz")] = nu
                row[key.replace("nu_t", "speed_mps")] = nu * lam / 2
            row["pulse_us"] = t_pulse * 1e6
            row["speed_per_unit_nut_mps"] = lam / (2 * t_pulse)
            summary.append(row)
            for j in range(0, cut["nu_t"].size, 5):  # every 0.05 in nu*T
                cuts.append(
                    {
                        "design_point": point,
                        "waveform": name,
                        "nu_t": cut["nu_t"][j],
                        "peak_loss_db": cut["peak_loss_db"][j],
                        "shift_cells": cut["shift_cells"][j],
                        "psl_db": cut["psl_db"][j],
                    }
                )
            print(
                f"{name:>14}: 1 dB at nuT {summ['nu_t_1db']:.3f}, 3 dB at {summ['nu_t_3db']:.3f}, "
                f"loss at nuT 0.5/1/2/4: {summ['loss_at_0p5_db']:.2f} {summ['loss_at_1_db']:.2f} "
                f"{summ['loss_at_2_db']:.2f} {summ['loss_at_4_db']:.2f} dB, "
                f"shift at 1: {summ['shift_at_1_cells']:+.2f} cells"
            )
    write_csv(OUT / "doppler_tolerance.csv", summary)
    write_csv(OUT / "doppler_cuts.csv", cuts)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
