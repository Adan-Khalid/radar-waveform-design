"""Low-bit ADC study: detection loss against bits, dither, ADC full scale and code length.

    python scripts/run_adc.py --study a1     # bits and dither, every waveform, BT 13 and 64
    python scripts/run_adc.py --study a1rep  # LFM BT 64 again with a new seed (replicate)
    python scripts/run_adc.py --study a2     # ADC full scale, mid-rise and mid-tread
    python scripts/run_adc.py --study a3     # time-bandwidth product 13 to 1024
    python scripts/run_adc.py --study pd     # wide Pd curves for the figure (LFM, BT 64)
    python scripts/run_adc.py --study phase  # one bit, 13 chips, fixed echo phase
    python scripts/run_adc.py --study a1 --quick --out <folder>   # tiny run for checks

Writes reports/adc/<study>_losses.csv, <study>_pd_curves.csv and <study>_pfa_check.csv.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from radarwave import waveforms as wf
from radarwave.config import REPORTS, load_config, write_csv
from radarwave.detection import Receiver, run_receivers

QUICK = {
    "pfa": 1e-3,
    "noise_statistics": 100_000,
    "trials_per_snr": 1500,
    "coarse_trials": 400,
    "bootstrap": 50,
}


def settings(cfg: dict, quick: bool) -> dict:
    c = dict(cfg["common"])
    out = {
        "pfa": c["pfa"],
        "pd_target": c["pd_target"],
        "noise_stats": c["noise_statistics"],
        "trials": c["trials_per_snr"],
        "coarse_trials": c["coarse_trials"],
        "bootstrap": c["bootstrap"],
        "fine_half_width": c["fine_half_width_db"],
        "fine_step": c["fine_step_db"],
    }
    if quick:
        out.update(
            pfa=QUICK["pfa"],
            noise_stats=QUICK["noise_statistics"],
            trials=QUICK["trials_per_snr"],
            coarse_trials=QUICK["coarse_trials"],
            bootstrap=QUICK["bootstrap"],
            fine_step=0.25,
        )
    return out


def run_set(s, receivers, opts, seed, tags: dict):
    loss, curves, pfa, _ = run_receivers(s, receivers, seed=seed, **opts)
    for rows in (loss, curves, pfa):
        for r in rows:
            r.update(tags)
    return loss, curves, pfa


def study_a1(cfg, opts, quick):
    a1 = cfg["a1_bits"]
    wcfg = load_config("waveforms")
    bits = [b for b in a1["bits"] if b is not None]
    dithers = a1["dither"]
    if quick:
        bits, dithers = [3, 1], ["none", "subtractive"]
    receivers = [Receiver()] + [
        Receiver(bits=b, kind=a1["quantizer"], dither=d) for b in bits for d in dithers
    ]
    out = ([], [], [])
    points = a1["design_points"] if not quick else ["long"]
    for point in points:
        dp = wcfg["design_points"][point]
        bt = dp["time_bandwidth"]
        specs = dp["waveforms"] if not quick else dp["waveforms"][:1]
        for spec in specs:
            s = wf.make_waveform(spec, bt, 1)
            tags = {
                "design_point": point,
                "time_bandwidth": bt,
                "waveform": wf.label(spec, bt),
                "papr_db": wf.papr_db(s),
            }
            t0 = time.time()
            res = run_set(s, receivers, opts, cfg["common"]["seed"], tags)
            for acc, rows in zip(out, res, strict=True):
                acc.extend(rows)
            _report(res[0], f"{tags['waveform']} ({time.time() - t0:.0f} s)")
    return out


def study_a1rep(cfg, opts, quick):
    """A1 for LFM BT 64 again with fresh random numbers (seed + 1): a full replicate,
    so the run-to-run spread, threshold included, can be compared with the intervals."""
    a1 = cfg["a1_bits"]
    bits = [b for b in a1["bits"] if b is not None] if not quick else [1]
    receivers = [Receiver()] + [
        Receiver(bits=b, kind=a1["quantizer"], dither=d) for b in bits for d in a1["dither"]
    ]
    s = wf.lfm(64, 1)
    seed = cfg["common"]["seed"] + 1
    tags = {"time_bandwidth": 64, "waveform": wf.label({"name": "lfm"}, 64), "seed": seed}
    res = run_set(s, receivers, opts, seed, tags)
    _report(res[0], "LFM (BT 64), replicate")
    return res


def study_a2(cfg, opts, quick):
    a2 = cfg["a2_full_scale"]
    fs = a2["full_scale_over_sigma"]
    grid_db = np.linspace(fs["start_db"], fs["stop_db"], fs["points"])
    bits, kinds, dithers = a2["bits"], a2["quantizer"], a2["dither"]
    if quick:
        grid_db, bits, kinds, dithers = np.array([0.0, 30.0]), [3], ["midrise"], ["none"]
    bt = a2["time_bandwidth"]
    s = wf.make_waveform(a2["waveform"], bt, 1)
    out = ([], [], [])
    for x_db in grid_db:
        a = float(10 ** (x_db / 20))
        group = [Receiver()]
        for b in bits:
            for kind in kinds:
                if kind == "midtread" and b < 2:
                    continue
                for d in dithers:
                    group.append(Receiver(bits=b, kind=kind, dither=d, full_scale=a))
        tags = {
            "time_bandwidth": bt,
            "waveform": wf.label(a2["waveform"], bt),
            "a_over_sigma_db": float(x_db),
        }
        t0 = time.time()
        res = run_set(s, group, opts, cfg["common"]["seed"], tags)
        for acc, rows in zip(out, res, strict=True):
            acc.extend(rows)
        print(f"A/sigma {x_db:6.1f} dB: {len(group)} receivers, {time.time() - t0:.0f} s")
    return out


def study_a3(cfg, opts, quick):
    a3 = cfg["a3_time_bandwidth"]
    bts = a3["time_bandwidth"] if not quick else [13, 64]
    bits = [b for b in a3["bits"] if b is not None] if not quick else [1]
    receivers = [Receiver()] + [Receiver(bits=b, kind=a3["quantizer"]) for b in bits]
    out = ([], [], [])
    for bt in bts:
        s = wf.make_waveform(a3["waveform"], bt, 1)
        tags = {"time_bandwidth": bt, "waveform": wf.label(a3["waveform"], bt)}
        t0 = time.time()
        res = run_set(s, receivers, opts, cfg["common"]["seed"], tags)
        for acc, rows in zip(out, res, strict=True):
            acc.extend(rows)
        _report(res[0], f"{tags['waveform']} ({time.time() - t0:.0f} s)")
    return out


def study_pd(cfg, opts, quick):
    """Full Pd-against-SNR curves (wide grid) for the figure: LFM BT 64, no dither."""
    from radarwave import theory
    from radarwave.detection import (
        TrialSet,
        check_pfa,
        detections,
        noise_statistics,
        thresholds_for,
        wilson,
    )

    bt = 64
    s = wf.lfm(bt, 1)
    seed = cfg["common"]["seed"]
    receivers = [Receiver()] + [Receiver(bits=b) for b in ([4, 3, 2, 1] if not quick else [1])]
    stats = noise_statistics(s, receivers, opts["noise_stats"], [seed, 1])
    thr = thresholds_for(stats, opts["pfa"])
    pfa_rows = check_pfa(s, receivers, thr, opts["noise_stats"] // 2, [seed, 2], opts["pfa"])
    grid = np.round(np.arange(4.0, 18.0 + 1e-9, 0.5 if not quick else 2.0), 6)
    trials = TrialSet(s.size, opts["trials"], [seed, 3])
    curves = []
    for rx in receivers:
        hits = detections(s, rx, thr[rx], grid, trials)
        lo, hi = wilson(hits.sum(axis=1), hits.shape[1])
        for j, snr in enumerate(grid):
            curves.append(
                {
                    "receiver": rx.label,
                    "bits": "" if rx.bits is None else rx.bits,
                    "snr_db": float(snr),
                    "pd": float(hits[j].mean()),
                    "pd_lo": float(lo[j]),
                    "pd_hi": float(hi[j]),
                    "trials": hits.shape[1],
                    "pd_marcum": float(theory.pd_marcum(snr, opts["pfa"])),
                    "waveform": wf.label({"name": "lfm"}, bt),
                }
            )
        print(f"{rx.label}: Pd curve over {grid[0]} to {grid[-1]} dB")
    return [], curves, pfa_rows


def study_phase(cfg, opts, quick):
    """One bit and a 13-chip code: Pd at a fixed echo phase instead of a random one.

    Shows why one-bit I/Q sampling of a short binary code loses so much (A1): a real code
    at phase 0 puts its whole echo into I. Uses the same receiver, threshold and decision
    code as every other study (`detections` with a fixed phase).
    """
    from radarwave.detection import (
        TrialSet,
        check_pfa,
        detections,
        noise_statistics,
        thresholds_for,
    )

    seed = cfg["common"]["seed"]
    phases = [0.0, 10.0, 22.5, 45.0] if not quick else [0.0, 45.0]
    snrs = [12.0, 16.0, 20.0, 25.0, 40.0] if not quick else [16.0, 25.0]
    rx = Receiver(bits=1)
    curves, pfa_rows = [], []
    for spec in ({"name": "barker", "length": 13}, {"name": "p4", "length": 13}):
        s = wf.make_waveform(spec, 13, 1)
        name = wf.label(spec, 13)
        thr = thresholds_for(noise_statistics(s, [rx], opts["noise_stats"], [seed, 1]), opts["pfa"])
        for row in check_pfa(s, [rx], thr, opts["noise_stats"] // 2, [seed, 2], opts["pfa"]):
            row["waveform"] = name
            pfa_rows.append(row)
        trials = TrialSet(s.size, opts["trials"], [seed, 3])
        for ph in phases:
            hits = detections(s, rx, thr[rx], snrs, trials, phase_deg=ph)
            for j, snr in enumerate(snrs):
                curves.append(
                    {
                        "waveform": name,
                        "receiver": rx.label,
                        "echo_phase_deg": ph,
                        "snr_db": snr,
                        "pd": float(hits[j].mean()),
                        "trials": hits.shape[1],
                    }
                )
            print(f"{name}, phase {ph:4.1f} deg: Pd {np.round(hits.mean(axis=1), 4).tolist()}")
    return [], curves, pfa_rows


def _report(rows, title):
    print(title)
    for r in rows:
        print(
            f"  {r['receiver']:>32}: req {r['required_snr_db']:6.2f} dB, "
            f"loss {r['loss_db']:5.2f} [{r['loss_lo']:5.2f}, {r['loss_hi']:5.2f}] "
            f"theory {r['theory_loss_db']:5.2f}"
        )


STUDIES = {
    "a1": study_a1,
    "a1rep": study_a1rep,
    "a2": study_a2,
    "a3": study_a3,
    "pd": study_pd,
    "phase": study_phase,
}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--study", choices=sorted(STUDIES), required=True)
    p.add_argument("--quick", action="store_true", help="tiny trial counts, for checks only")
    p.add_argument("--out", type=Path, default=None, help="output folder (default reports/adc)")
    args = p.parse_args(argv)
    cfg = load_config("adc")
    opts = settings(cfg, args.quick)
    out_dir = args.out or REPORTS / "adc"
    t0 = time.time()
    loss, curves, pfa = STUDIES[args.study](cfg, opts, args.quick)
    for name, rows in (("losses", loss), ("pd_curves", curves), ("pfa_check", pfa)):
        if rows:  # some studies have no loss table or no false-alarm check
            write_csv(out_dir / f"{args.study}_{name}.csv", rows)
    meta = [{"setting": k, "value": v} for k, v in opts.items()]
    meta += [
        {"setting": "seed", "value": cfg["common"]["seed"]},
        {"setting": "quick", "value": args.quick},
    ]
    meta += [{"setting": "runtime_s", "value": round(time.time() - t0, 1)}]
    write_csv(out_dir / f"{args.study}_settings.csv", meta)
    print(f"wrote {out_dir} in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
