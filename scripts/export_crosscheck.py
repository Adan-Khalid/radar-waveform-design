"""Reference values for the MATLAB cross-check (small CSV files in reports/crosscheck/).

MATLAB's random generator differs from NumPy's, so the OFDM symbols are exported and
MATLAB builds its OFDM pulses from them. Waveform samples and filter taps are exported so
MATLAB can check its own versions against them sample by sample.
"""

from __future__ import annotations

from radarwave import waveforms as wf
from radarwave.config import REPORTS, load_config, write_csv
from radarwave.sidelobes import ls_mismatched_filter

OUT = REPORTS / "crosscheck"


def main():
    cfg = load_config("waveforms")
    rows = []
    for dp in cfg["design_points"].values():
        for spec in dp["waveforms"]:
            if spec["name"] != "ofdm":
                continue
            sym = wf.qpsk_symbols(spec["subcarriers"], spec["seed"])
            for k, x in enumerate(sym):
                rows.append(
                    {"subcarriers": spec["subcarriers"], "index": k, "re": x.real, "im": x.imag}
                )
    write_csv(OUT / "ofdm_symbols.csv", rows, digits=17)

    samples = []
    cases = [
        ("lfm_bt13_spc8", wf.lfm(13, 8)),
        ("p4_13_spc8", wf.phase_code(wf.p4_phases(13), 8)),
        ("frank_64_spc1", wf.phase_code(wf.frank_phases(64), 1)),
        ("ofdm_13_spc2", wf.ofdm(wf.qpsk_symbols(13, 2026), 2)),
    ]
    for name, s in cases:
        for n, x in enumerate(s):
            samples.append({"signal": name, "n": n, "re": x.real, "im": x.imag})
    write_csv(OUT / "waveform_samples.csv", samples, digits=17)

    code = wf.phase_code(wf.barker_phases(13))
    ref = ls_mismatched_filter(code, 39, 1e-9)
    taps = [{"n": n, "re": x.real, "im": x.imag} for n, x in enumerate(ref)]
    write_csv(OUT / "ls_filter_barker13_x3.csv", taps, digits=17)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
