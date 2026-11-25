"""Sidelobe control and its SNR cost: windows on LFM, mismatched filters on phase codes.

Writes reports/sidelobes/lfm_windows.csv, mismatched_filters.csv and mismatched_doppler.csv.
"""

from __future__ import annotations

import numpy as np

from radarwave import theory
from radarwave import waveforms as wf
from radarwave.ambiguity import doppler_cut_summary
from radarwave.config import REPORTS, load_config, write_csv
from radarwave.metrics import pulse_metrics
from radarwave.sidelobes import (
    ls_mismatched_filter,
    minpsl_filter_real,
    window,
    window_label,
    windowed_reference,
)

OUT = REPORTS / "sidelobes"


def main():
    cfg = load_config("sidelobes")
    spc = load_config("waveforms")["analysis"]["samples_per_cell"]

    rows = []
    for bt in cfg["lfm_windows"]["time_bandwidth"]:
        s = wf.lfm(bt, spc)
        for spec in cfg["lfm_windows"]["windows"]:
            kw = {k: v for k, v in spec.items() if k != "name"}
            w = window(spec["name"], s.size, **kw)
            m = pulse_metrics(s, windowed_reference(s, w), spc)
            rows.append(
                {
                    "time_bandwidth": bt,
                    "window": window_label(spec),
                    **m,
                    "window_loss_formula_db": theory.window_loss_db(w),
                }
            )
            print(
                f"LFM BT {bt:>3} {window_label(spec):>20}: PSL {m['psl_db']:7.2f} "
                f"ISL {m['isl_db']:7.2f} width {m['width_3db_cells']:.3f} "
                f"loss {m['mismatch_loss_db']:.3f} dB"
            )
    write_csv(OUT / "lfm_windows.csv", rows)

    mm = cfg["mismatched_filter"]
    rows, dop = [], []
    nu = np.round(np.arange(0, 1.0 + 1e-9, 0.01), 10)
    for spec in mm["codes"]:
        n = spec["length"]
        code = wf.make_waveform(spec, n, 1)  # one sample per chip
        name = wf.label(spec, n)
        filters = [("matched", 1, code)]
        for f in mm["length_factors"]:
            ref = ls_mismatched_filter(code, f * n, mm["diagonal_loading"])
            filters.append(("least squares", f, ref))
        if spec["name"] == "barker":
            filters.append(("minimum PSL (LP)", 3, minpsl_filter_real(code, 3 * n)))
        for kind, f, ref in filters:
            m = pulse_metrics(code, ref, 1)
            rows.append(
                {
                    "code": name,
                    "filter": kind,
                    "length_factor": f,
                    "filter_length": ref.size,
                    **m,
                }
            )
            print(
                f"{name:>10} {kind:>17} x{f}: PSL {m['psl_db']:7.2f} ISL {m['isl_db']:7.2f} "
                f"loss {m['mismatch_loss_db']:.3f} dB"
            )
            if kind == "matched" or (kind == "least squares" and f == 3):
                cut = doppler_cut_summary(code, ref, nu, 1)
                for j in range(nu.size):
                    dop.append(
                        {
                            "code": name,
                            "filter": kind,
                            "nu_t": nu[j],
                            "peak_loss_db": cut["peak_loss_db"][j],
                            "psl_db": cut["psl_db"][j],
                        }
                    )
    write_csv(OUT / "mismatched_filters.csv", rows)
    write_csv(OUT / "mismatched_doppler.csv", dop)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
