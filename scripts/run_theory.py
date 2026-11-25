"""Theory tables: required SNR, window losses and quantizer losses (Bussgang efficiency).

Writes reports/theory/*.csv. No randomness.
"""

from __future__ import annotations

import numpy as np

from radarwave import theory
from radarwave.adc import GAUSSIAN_OPTIMAL_STEP, gaussian_optimal_full_scale
from radarwave.config import REPORTS, load_config, write_csv
from radarwave.sidelobes import window, window_label

OUT = REPORTS / "theory"


def main():
    adc = load_config("adc")
    pfa = adc["common"]["pfa"]
    pd = adc["common"]["pd_target"]

    rows = []
    for p in (1e-3, 1e-4, 1e-6):
        for d in (0.5, 0.9):
            rows.append(
                {
                    "pfa": p,
                    "pd": d,
                    "required_snr_marcum_db": theory.required_snr_db(d, p),
                    "required_snr_albersheim_db": theory.albersheim_db(d, p),
                }
            )
    write_csv(OUT / "required_snr.csv", rows)

    wins = load_config("sidelobes")["lfm_windows"]["windows"]
    rows = []
    for spec in wins:
        kw = {k: v for k, v in spec.items() if k != "name"}
        w = window(spec["name"], 100_001, **kw)
        rows.append({"window": window_label(spec), "loss_db": theory.window_loss_db(w)})
    write_csv(OUT / "window_loss.csv", rows)

    # Quantizer loss at Max's Gaussian-optimal full scale, every kind and dither.
    rows = []
    for bits in (1, 2, 3, 4):
        a = gaussian_optimal_full_scale(bits)
        for kind in ("midrise", "midtread"):
            if kind == "midtread" and bits == 1:
                continue
            for dither in ("none", "nonsubtractive", "subtractive"):
                rows.append(
                    {
                        "bits": bits,
                        "kind": kind,
                        "dither": dither,
                        "step_over_sigma": GAUSSIAN_OPTIMAL_STEP[bits],
                        "a_over_sigma": a,
                        "efficiency": theory.quantizer_efficiency(bits, a, kind, dither),
                        "loss_db": theory.quantizer_loss_db(bits, a, kind, dither),
                    }
                )
    write_csv(OUT / "quantizer_loss_optimal.csv", rows)

    # Loss against full scale (A / sigma per component), for the full-scale study figure.
    a_db = np.linspace(-15, 50, 261)
    rows = []
    for bits in (1, 2, 3, 4):
        for kind in ("midrise", "midtread"):
            if kind == "midtread" and bits == 1:
                continue
            for dither in ("none", "nonsubtractive", "subtractive"):
                for x in a_db:
                    a = 10 ** (x / 20)
                    rows.append(
                        {
                            "bits": bits,
                            "kind": kind,
                            "dither": dither,
                            "a_over_sigma_db": x,
                            "loss_db": theory.quantizer_loss_db(bits, a, kind, dither),
                        }
                    )
    write_csv(OUT / "quantizer_loss_vs_full_scale.csv", rows)
    print(f"required SNR (Pd {pd}, Pfa {pfa}): {theory.required_snr_db(pd, pfa):.3f} dB")
    for bits in (1, 2, 3, 4):
        a = gaussian_optimal_full_scale(bits)
        print(f"{bits} bit at Max's step: {theory.quantizer_loss_db(bits, a):.4f} dB")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
