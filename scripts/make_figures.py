"""Every figure in reports/figures/, from the CSV files in reports/ (and the deterministic
waveforms, which are cheap to rebuild). Run the study scripts first.

    python scripts/make_figures.py            # all figures whose inputs exist
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from radarwave import theory
from radarwave import waveforms as wf
from radarwave.ambiguity import ambiguity
from radarwave.config import REPORTS, load_config, read_csv, use_style
from radarwave.matched import correlate
from radarwave.sidelobes import ls_mismatched_filter, minpsl_filter_real, window

plt = use_style()
FIG = REPORTS / "figures"
FIG.mkdir(parents=True, exist_ok=True)

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
# Waveform families keep one colour in every figure (categorical slots 1 to 4).
FAMILY = {
    "lfm": "#2a78d6",
    "barker": "#eb6834",
    "frank": "#eb6834",
    "p4": "#1baf7a",
    "ofdm": "#eda100",
}
# Bit depth is ordered, so it uses one blue ramp: more bits, lighter.
BITS = {1: "#0d366b", 2: "#1c5cab", 3: "#3987e5", 4: "#86b6ef"}
DITHER_STYLE = {
    "none": dict(color=INK, marker="o", ls="-", label="no dither"),
    "nonsubtractive": dict(color=MUTED, marker="s", ls="--", label="dither (not removed)"),
    "subtractive": dict(color=INK2, marker="^", ls=":", label="dither (removed after ADC)"),
}


def family(name: str) -> str:
    return name.split()[0].lower()


def save(fig, name):
    path = FIG / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def db(x, floor=-80.0):
    return np.maximum(20 * np.log10(np.maximum(np.abs(x), 1e-12)), floor)


def _place_labels(ax, ends, min_gap):
    """Direct labels at line ends, pushed apart vertically so they never overlap."""
    ends = sorted(ends, key=lambda e: e[0], reverse=True)
    for i in range(1, len(ends)):
        if ends[i - 1][0] - ends[i][0] < min_gap:
            ends[i][0] = ends[i - 1][0] - min_gap
    for y, name, x in ends:
        ax.annotate(
            name,
            (x, y),
            xytext=(4, 0),
            textcoords="offset points",
            fontsize=7.5,
            color=INK2,
            va="center",
        )


# ---------------------------------------------------------------- pulse compression


def fig_compressed_pulses():
    cfg = load_config("waveforms")
    spc = cfg["analysis"]["samples_per_cell"]
    metrics = {
        (r["design_point"], r["waveform"]): r for r in read_csv(REPORTS / "waveforms/metrics.csv")
    }
    points = list(cfg["design_points"].items())
    fig, axes = plt.subplots(2, 4, figsize=(11, 4.6), sharey=True)
    for i, (point, dp) in enumerate(points):
        bt = dp["time_bandwidth"]
        for j, spec in enumerate(dp["waveforms"]):
            ax = axes[i, j]
            s = wf.make_waveform(spec, bt, spc)
            y, lag = correlate(s, s)
            name = wf.label(spec, bt)
            ax.plot(lag / spc, db(y), color=FAMILY[spec["name"]], lw=1.0)
            m = metrics[(point, name)]
            ax.set_title(f"{name}", fontsize=9)
            ax.text(
                0.03,
                0.95,
                f"PSL {float(m['psl_db']):.1f} dB\nISL {float(m['isl_db']):.1f} dB",
                transform=ax.transAxes,
                va="top",
                fontsize=7.5,
                color=INK2,
            )
            ax.set_xlim(-bt, bt)
            ax.set_ylim(-60, 3)
            if j == 0:
                ax.set_ylabel(f"BT = {bt}\nmatched-filter output (dB)")
            if i == 1:
                ax.set_xlabel("delay (range cells of 1/B)")
    fig.suptitle(
        "Compressed pulses at equal bandwidth and duration", x=0.01, ha="left", fontweight="bold"
    )
    fig.tight_layout()
    save(fig, "compressed_pulses")


def fig_ambiguity():
    cfg = load_config("waveforms")
    spc = 4
    dp = cfg["design_points"]["long"]
    bt = dp["time_bandwidth"]
    nu = np.linspace(-4, 4, 241)
    max_lag = 8
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.5), sharey=True)
    fig.subplots_adjust(top=0.82, bottom=0.16, left=0.07, right=0.9, wspace=0.12)
    for ax, spec in zip(axes, dp["waveforms"], strict=True):
        s = wf.make_waveform(spec, bt, spc)
        lags = np.arange(-max_lag * spc, max_lag * spc + 1)
        chi, _, _ = ambiguity(s, nu_t=nu, lags=lags)
        im = ax.imshow(
            db(chi).T,
            extent=[lags[0] / spc, lags[-1] / spc, nu[0], nu[-1]],
            aspect="auto",
            cmap="Blues",
            vmin=-30,
            vmax=0,
            origin="lower",
        )
        ax.set_title(wf.label(spec, bt), fontsize=9)
        ax.set_xlabel("delay (range cells)")
        ax.grid(False)
    axes[0].set_ylabel("Doppler x pulse length, νT")
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.01)
    cb.set_label("|χ| (dB)")
    fig.suptitle(
        "Ambiguity functions, BT = 64 (darker = stronger response)",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    save(fig, "ambiguity_functions")


def fig_doppler_tolerance():
    rows = read_csv(REPORTS / "ambiguity/doppler_cuts.csv")
    data = defaultdict(list)
    for r in rows:
        data[(r["design_point"], r["waveform"])].append(
            (float(r["nu_t"]), float(r["peak_loss_db"]))
        )
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    nu = np.linspace(0, 4, 401)
    for ax, point, bt in zip(axes, ("short", "long"), (13, 64), strict=True):
        ax.plot(nu, db(np.sinc(nu)), color=MUTED, lw=1, ls="--")
        ax.text(2.15, -15.2, "|sinc(νT)|", color=MUTED, fontsize=7.5)
        ends = []
        for (p, name), pts in data.items():
            if p != point:
                continue
            x, y = np.array(pts).T
            ax.plot(x, y, color=FAMILY[family(name)], lw=1.5)
            ends.append([y[-1], name, x[-1]])
        _place_labels(ax, ends, min_gap=0.9)
        ax.set_xlim(0, 4.9)
        ax.set_ylim(-20, 1)
        ax.set_title(f"BT = {bt}", fontsize=9)
        ax.set_xlabel("Doppler x pulse length, νT")
    axes[0].set_ylabel("peak of the compressed pulse (dB)")
    fig.suptitle(
        "Doppler tolerance: how much of the peak survives a Doppler shift",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout()
    save(fig, "doppler_tolerance")


def fig_sidelobe_control():
    spc = load_config("waveforms")["analysis"]["samples_per_cell"]
    win = read_csv(REPORTS / "sidelobes/lfm_windows.csv")
    mm = read_csv(REPORTS / "sidelobes/mismatched_filters.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), gridspec_kw={"width_ratios": [1, 1, 1.15]})

    # (a) LFM BT 64 with windows
    ax = axes[0]
    s = wf.lfm(64, spc)
    styles = [
        ("rectangular", {}, INK2, "-"),
        ("taylor", {"sidelobe_db": 35, "nbar": 5}, FAMILY["lfm"], "-"),
        ("hamming", {}, FAMILY["p4"], "--"),
    ]
    for name, kw, col, ls in styles:
        y, lag = correlate(s, s * window(name, s.size, **kw))
        y = y / np.max(np.abs(y))
        ax.plot(
            lag / spc,
            db(y),
            color=col,
            ls=ls,
            lw=1.2,
            label=name.capitalize() if name != "taylor" else "Taylor 35 dB",
        )
    ax.set_xlim(-12, 12)
    ax.set_ylim(-60, 3)
    ax.set_xlabel("delay (range cells)")
    ax.set_ylabel("output (dB, own peak = 0)")
    ax.set_title("LFM, BT 64: receive windows", fontsize=9)
    ax.legend(loc="upper right")

    # (b) Barker 13: matched, least squares x3, minimum PSL x3
    ax = axes[1]
    code = wf.phase_code(wf.barker_phases(13))
    for ref, col, ls, lab in [
        (code, INK2, "-", "matched"),
        (ls_mismatched_filter(code, 39), FAMILY["barker"], "-", "least squares, 39 taps"),
        (minpsl_filter_real(code, 39), FAMILY["lfm"], "--", "minimum PSL, 39 taps"),
    ]:
        y, lag = correlate(code, ref)
        y = y / np.max(np.abs(y))
        ax.plot(
            lag - lag[np.argmax(np.abs(y))],
            db(y),
            color=col,
            ls=ls,
            lw=1.2,
            marker=".",
            ms=3,
            label=lab,
        )
    ax.set_xlim(-30, 30)
    ax.set_ylim(-70, 3)
    ax.set_xlabel("delay (chips)")
    ax.set_title("Barker 13: mismatched filters", fontsize=9)
    ax.legend(loc="upper right")

    # (c) trade-off: PSL against SNR cost, every option
    ax = axes[2]
    for r in win:
        if int(r["time_bandwidth"]) != 64:
            continue
        x, y = float(r["mismatch_loss_db"]), float(r["psl_db"])
        ax.scatter(x, y, color=FAMILY["lfm"], marker="o", s=28, zorder=3)
        ax.annotate(
            r["window"].split(",")[0],
            (x, y),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=7,
            color=INK2,
        )
    for code_name, mk in (("Barker 13", "s"), ("Frank 64", "D"), ("P4 64", "^")):
        pts = [
            (float(r["mismatch_loss_db"]), float(r["psl_db"]), r["filter"], r["length_factor"])
            for r in mm
            if r["code"] == code_name
        ]
        col = FAMILY[family(code_name)]
        ls_pts = sorted(
            (p for p in pts if p[2] in ("matched", "least squares")), key=lambda p: p[0]
        )
        ax.plot(
            [p[0] for p in ls_pts],
            [p[1] for p in ls_pts],
            color=col,
            lw=1,
            marker=mk,
            ms=5,
            label=f"{code_name}: matched, then LS x1 to x4",
        )
        for p in pts:
            if p[2].startswith("minimum"):
                ax.scatter(
                    p[0],
                    p[1],
                    color=col,
                    marker="*",
                    s=70,
                    zorder=4,
                    label=f"{code_name}: minimum PSL x3",
                )
    ax.set_xlabel("SNR loss against the matched filter (dB)")
    ax.set_ylabel("peak sidelobe level (dB)")
    ax.set_title("Sidelobe level against its SNR cost (LFM windows: blue dots)", fontsize=9)
    ax.legend(loc="upper right", fontsize=6.5)
    ax.set_xlim(-0.1, 2.6)
    fig.suptitle("Sidelobe control and what it costs", x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    save(fig, "sidelobe_control")


def fig_mismatched_doppler():
    rows = read_csv(REPORTS / "sidelobes/mismatched_doppler.csv")
    data = defaultdict(list)
    for r in rows:
        data[(r["code"], r["filter"])].append((float(r["nu_t"]), float(r["psl_db"])))
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for ax, codes in zip(axes, (("Barker 13", "P4 13"), ("Frank 64", "P4 64")), strict=True):
        for (code, filt), pts in data.items():
            if code not in codes:
                continue
            x, y = np.array(pts).T
            ls = "-" if filt == "matched" else "--"
            tag = "matched filter" if filt == "matched" else "least squares, 3N taps"
            ax.plot(x, y, color=FAMILY[family(code)], ls=ls, lw=1.4, label=f"{code}, {tag}")
        ax.set_xlabel("Doppler x pulse length, νT")
        ax.set_xlim(0, 0.5)
        ax.legend(fontsize=7, loc="lower right")
    axes[0].set_ylabel("PSL relative to the peak (dB)")
    fig.suptitle(
        "Mismatched filters lose their low sidelobes as Doppler grows",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout()
    save(fig, "mismatched_doppler")


# ---------------------------------------------------------------- low-bit ADC


def _f(x):
    return float(x) if x not in ("", None) else float("nan")


def fig_adc_bits():
    rows = read_csv(REPORTS / "adc/a1_losses.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"width_ratios": [1.15, 1]})
    # (a) LFM BT 64: loss against bits for each dither setting, with theory
    ax = axes[0]
    sel = [r for r in rows if r["waveform"] == "LFM (BT 64)" and r["bits"] != ""]
    offsets = {"none": -0.12, "nonsubtractive": 0.0, "subtractive": 0.12}
    for dither, st in DITHER_STYLE.items():
        pts = sorted((int(r["bits"]), r) for r in sel if r["dither"] == dither)
        b = np.array([p[0] for p in pts]) + offsets[dither]
        loss = np.array([_f(p[1]["loss_db"]) for p in pts])
        lo = np.array([_f(p[1]["loss_lo"]) for p in pts])
        hi = np.array([_f(p[1]["loss_hi"]) for p in pts])
        th = np.array([_f(p[1]["theory_loss_db"]) for p in pts])
        ax.plot(b, th, color=st["color"], ls=st["ls"], lw=1)
        ax.errorbar(
            b,
            loss,
            yerr=[loss - lo, hi - loss],
            fmt=st["marker"],
            color=st["color"],
            ms=6,
            capsize=2.5,
            lw=1,
            label=st["label"],
        )
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xlabel("ADC bits (I and Q each)")
    ax.set_ylabel("detection loss against unquantized (dB)")
    ax.set_title(
        "LFM, BT 64: markers = Monte Carlo with 95 % interval, lines = theory", fontsize=8.5
    )
    ax.legend(loc="upper right")
    ax.set_ylim(bottom=0)
    # (b) no dither: every waveform at BT 13 (open) and BT 64 (filled)
    ax = axes[1]
    bb = np.arange(1, 5)
    th = [theory.quantizer_loss_db(b, _gauss_a(b)) for b in bb]
    ax.plot(bb, th, color=MUTED, ls="--", lw=1)
    ax.text(3.05, th[2] + 0.35, "theory\n(low SNR per sample)", color=MUTED, fontsize=7)
    names = sorted({r["waveform"] for r in rows})
    for name in names:
        pts = sorted(
            (int(r["bits"]), _f(r["loss_db"]))
            for r in rows
            if r["waveform"] == name and r["dither"] == "none" and r["bits"] != ""
        )
        bt = int(next(r["time_bandwidth"] for r in rows if r["waveform"] == name))
        x = np.array([p[0] for p in pts]) + (0.06 if bt == 64 else -0.06)
        col = FAMILY[family(name)]
        ax.plot(
            x,
            [p[1] for p in pts],
            ls="none",
            marker="o",
            ms=5.5,
            mfc=col if bt == 64 else "none",
            mec=col,
            label=name,
        )
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xlabel("ADC bits (no dither)")
    ax.set_title("Every waveform: BT 13 (open) and BT 64 (filled)", fontsize=8.5)
    ax.legend(fontsize=6.5, ncol=2, loc="upper right")
    ax.set_ylim(bottom=0)
    fig.suptitle(
        "What each ADC bit costs in detection (Pd 0.9, Pfa 1e-4)",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout()
    save(fig, "adc_loss_vs_bits")


def _gauss_a(bits):
    from radarwave.adc import gaussian_optimal_full_scale

    return gaussian_optimal_full_scale(bits)


def fig_adc_pd_curves():
    path = REPORTS / "adc/pd_pd_curves.csv"
    if not path.exists():
        return
    rows = read_csv(path)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    snr = np.linspace(4, 18, 300)
    ax.plot(
        snr,
        theory.pd_marcum(snr, 1e-4),
        color=MUTED,
        ls="--",
        lw=1.2,
        label="theory, unquantized (Marcum Q)",
    )
    groups = defaultdict(list)
    for r in rows:
        groups[r["receiver"]].append(r)
    for label, pts in groups.items():
        x = np.array([_f(p["snr_db"]) for p in pts])
        y = np.array([_f(p["pd"]) for p in pts])
        lo = np.array([_f(p["pd_lo"]) for p in pts])
        hi = np.array([_f(p["pd_hi"]) for p in pts])
        bits = pts[0]["bits"]
        col = INK if bits == "" else BITS[int(bits)]
        ax.fill_between(x, lo, hi, color=col, alpha=0.18, lw=0)
        ax.plot(x, y, color=col, marker="o", ms=3, lw=1.2, label=label)
    ax.axhline(0.9, color=MUTED, lw=0.6)
    ax.set_xlabel("SNR after the matched filter, unquantized (dB)")
    ax.set_ylabel("probability of detection")
    ax.set_title("Pd at Pfa 1e-4, LFM BT 64, no dither (bands: 95 % intervals)", fontsize=9)
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    save(fig, "adc_pd_curves")


def fig_adc_full_scale():
    path = REPORTS / "adc/a2_losses.csv"
    if not path.exists():
        return
    rows = [r for r in read_csv(path) if r["bits"] != ""]
    th = read_csv(REPORTS / "theory/quantizer_loss_vs_full_scale.csv")
    panels = [
        ("midrise", "none"),
        ("midrise", "subtractive"),
        ("midrise", "nonsubtractive"),
        ("midtread", "none"),
    ]
    titles = {
        ("midrise", "none"): "mid-rise, no dither",
        ("midrise", "subtractive"): "mid-rise, dither removed after the ADC",
        ("midrise", "nonsubtractive"): "mid-rise, dither not removed",
        ("midtread", "none"): "mid-tread (has a zero level), no dither",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.4), sharex=True, sharey=True)
    for ax, (kind, dither) in zip(axes.ravel(), panels, strict=True):
        for bits in (1, 2, 3, 4):
            if kind == "midtread" and bits == 1:
                continue
            col = BITS[bits]
            t = [
                (float(r["a_over_sigma_db"]), _f(r["loss_db"]))
                for r in th
                if int(r["bits"]) == bits and r["kind"] == kind and r["dither"] == dither
            ]
            tx, ty = np.array(t).T
            ax.plot(tx, np.minimum(ty, 30), color=col, lw=1.2)
            pts = sorted(
                (_f(r["a_over_sigma_db"]), _f(r["loss_db"]), _f(r["loss_lo"]), _f(r["loss_hi"]))
                for r in rows
                if int(r["bits"]) == bits and r["kind"] == kind and r["dither"] == dither
            )
            if pts:
                x, y, lo, hi = np.array(pts).T
                ok = np.isfinite(y)
                ax.errorbar(
                    x[ok],
                    y[ok],
                    yerr=[y[ok] - lo[ok], hi[ok] - y[ok]],
                    fmt="o",
                    ms=3.5,
                    color=col,
                    lw=0.8,
                    capsize=1.5,
                    label=f"{bits} bit",
                )
        ax.set_title(titles[(kind, dither)], fontsize=9)
        ax.set_ylim(-0.2, 12)
        ax.set_xlim(-13, 47)
        ax.legend(loc="upper left", fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("ADC full scale over noise, A/σ per component (dB)")
    for ax in axes[:, 0]:
        ax.set_ylabel("detection loss (dB)")
    sec = axes[0, 0].secondary_xaxis("top", functions=(lambda x: x - 3.01, lambda x: x + 3.01))
    sec.set_xlabel(
        "strong echo setting the gain: SNR per sample, constant envelope (dB)", fontsize=7.5
    )
    fig.suptitle(
        "Loss against ADC full scale (LFM BT 64): dots = Monte Carlo with 95 % interval, "
        "lines = theory",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout()
    save(fig, "adc_full_scale")


def fig_adc_time_bandwidth():
    path = REPORTS / "adc/a3_losses.csv"
    if not path.exists():
        return
    rows = [r for r in read_csv(path) if r["bits"] != ""]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for bits in (1, 2, 3, 4):
        pts = sorted(
            (int(r["time_bandwidth"]), _f(r["loss_db"]), _f(r["loss_lo"]), _f(r["loss_hi"]))
            for r in rows
            if int(r["bits"]) == bits
        )
        if not pts:
            continue
        x, y, lo, hi = np.array(pts).T
        col = BITS[bits]
        ax.errorbar(
            x,
            y,
            yerr=[y - lo, hi - y],
            fmt="o-",
            color=col,
            ms=4,
            lw=1.2,
            capsize=2,
            label=f"{bits} bit",
        )
        ax.axhline(theory.quantizer_loss_db(bits, _gauss_a(bits)), color=col, ls="--", lw=0.8)
    ax.set_xscale("log", base=2)
    ax.set_xticks([13, 32, 64, 128, 256, 1024])
    ax.set_xticklabels(["13", "32", "64", "128", "256", "1024"])
    ax.set_xlabel("time-bandwidth product BT (P4 code, N = BT chips)")
    ax.set_ylabel("detection loss against unquantized (dB)")
    ax.set_title("Short codes lose more than the low-SNR theory (dashed)", fontsize=9)
    ax.legend(loc="upper right")
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    save(fig, "adc_time_bandwidth")


def main():
    fig_compressed_pulses()
    fig_ambiguity()
    fig_doppler_tolerance()
    fig_sidelobe_control()
    fig_mismatched_doppler()
    if (REPORTS / "adc/a1_losses.csv").exists():
        fig_adc_bits()
    fig_adc_pd_curves()
    fig_adc_full_scale()
    fig_adc_time_bandwidth()


if __name__ == "__main__":
    main()
