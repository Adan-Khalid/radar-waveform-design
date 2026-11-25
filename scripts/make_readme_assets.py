"""README assets: the headline numbers table, the hero image and the social preview.

Reads only files in reports/ (run the study scripts first) and writes:
- reports/summary/headline.csv   every number quoted in the README's "Results at a glance"
- assets/hero.png                the picture at the top of the README
- assets/social-preview.png      1280 x 640 image for the repository settings
"""

from __future__ import annotations

import numpy as np

from radarwave.config import REPO_ROOT, REPORTS, read_csv, use_style, write_csv

plt = use_style()
ASSETS = REPO_ROOT / "assets"
TITLE = (
    "Radar Waveform Design and Matched Filtering: "
    "LFM, Phase-Coded and OFDM Pulses with Low-Bit ADC Analysis"
)
INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#898781", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def _f(x):
    return float(x) if x not in ("", None) else float("nan")


def _row(rows, **match):
    hits = [r for r in rows if all(r.get(k) == v for k, v in match.items())]
    if len(hits) != 1:
        raise LookupError(f"{len(hits)} rows match {match}")
    return hits[0]


def headline():
    a1 = read_csv(REPORTS / "adc/a1_losses.csv")
    a2 = read_csv(REPORTS / "adc/a2_losses.csv")
    a3 = read_csv(REPORTS / "adc/a3_losses.csv")
    mm = read_csv(REPORTS / "sidelobes/mismatched_filters.csv")
    dt = read_csv(REPORTS / "ambiguity/doppler_tolerance.csv")
    th = read_csv(REPORTS / "theory/required_snr.csv")
    out = []

    def add(key, value, lo="", hi="", theory="", source="", where=""):
        out.append(
            {
                "key": key,
                "value": value,
                "lo": lo,
                "hi": hi,
                "theory": theory,
                "source": source,
                "row": where,
            }
        )

    lfm = "LFM (BT 64)"
    unq = _row(a1, waveform=lfm, receiver="unquantized")
    marcum = _row(th, pfa="0.0001", pd="0.9")
    add(
        "required_snr_unquantized_db",
        _f(unq["required_snr_db"]),
        _f(unq["required_snr_lo"]),
        _f(unq["required_snr_hi"]),
        _f(marcum["required_snr_marcum_db"]),
        "reports/adc/a1_losses.csv",
        f"{lfm}, unquantized",
    )
    for bits in ("4", "3", "2", "1"):
        r = _row(a1, waveform=lfm, receiver=f"{bits} bit")
        add(
            f"loss_{bits}bit_bt64_db",
            _f(r["loss_db"]),
            _f(r["loss_lo"]),
            _f(r["loss_hi"]),
            _f(r["theory_loss_db"]),
            "reports/adc/a1_losses.csv",
            f"{lfm}, {bits} bit",
        )
    for bt, bits in (("13", "1"), ("1024", "1"), ("1024", "3")):
        r = _row(a3, time_bandwidth=bt, receiver=f"{bits} bit")
        add(
            f"loss_{bits}bit_p4_bt{bt}_db",
            _f(r["loss_db"]),
            _f(r["loss_lo"]),
            _f(r["loss_hi"]),
            _f(r["theory_loss_db"]),
            "reports/adc/a3_losses.csv",
            f"P4 {bt}, {bits} bit",
        )
    for label, kind, dither in (
        ("midrise_nodither", "midrise", "none"),
        ("midrise_dither_removed", "midrise", "subtractive"),
        ("midtread_nodither", "midtread", "none"),
    ):
        r = _row(a2, a_over_sigma_db="30", bits="3", kind=kind, dither=dither)
        add(
            f"loss_3bit_fullscale30db_{label}",
            _f(r["loss_db"]),
            _f(r["loss_lo"]),
            _f(r["loss_hi"]),
            _f(r["theory_loss_db"]),
            "reports/adc/a2_losses.csv",
            f"A/sigma 30 dB, 3 bit, {kind}, {dither}",
        )
    for filt, factor in (("matched", "1"), ("least squares", "3")):
        r = _row(mm, code="Barker 13", filter=filt, length_factor=factor)
        add(
            f"barker13_{filt.replace(' ', '_')}_psl_db",
            _f(r["psl_db"]),
            source="reports/sidelobes/mismatched_filters.csv",
            where=f"Barker 13, {filt} x{factor}",
        )
        add(
            f"barker13_{filt.replace(' ', '_')}_loss_db",
            _f(r["mismatch_loss_db"]),
            source="reports/sidelobes/mismatched_filters.csv",
            where=f"Barker 13, {filt} x{factor}",
        )
    for wf in ("Frank 64", "P4 64", "LFM (BT 64)"):
        r = _row(dt, waveform=wf)
        add(
            f"doppler_3db_nut_{wf.split()[0].lower()}64",
            _f(r["nu_t_3db"]),
            source="reports/ambiguity/doppler_tolerance.csv",
            where=wf,
        )
    write_csv(REPORTS / "summary/headline.csv", out)
    print(f"wrote {REPORTS / 'summary/headline.csv'}")
    return {r["key"]: r for r in out}


def _bits_panel(ax, label_size=8):
    """Loss against bits for the P4 code at three lengths, with the low-SNR theory."""
    a3 = read_csv(REPORTS / "adc/a3_losses.csv")
    bits = np.arange(1, 5)
    theory = [
        _f(_row(a3, time_bandwidth="1024", receiver=f"{b} bit")["theory_loss_db"]) for b in bits
    ]
    ax.plot(bits, theory, color=INK, ls=":", lw=1.6, zorder=4, label="theory (low SNR per sample)")
    styles = {"13": (ORANGE, "s"), "64": (BLUE, "o"), "1024": (AQUA, "D")}
    for bt, (col, mk) in styles.items():
        rows = [_row(a3, time_bandwidth=bt, receiver=f"{b} bit") for b in bits]
        y = np.array([_f(r["loss_db"]) for r in rows])
        lo = np.array([_f(r["loss_lo"]) for r in rows])
        hi = np.array([_f(r["loss_hi"]) for r in rows])
        ax.errorbar(
            bits,
            y,
            yerr=[y - lo, hi - y],
            color=col,
            marker=mk,
            ms=6,
            lw=1.6,
            capsize=2,
            zorder=3,
            label=f"P4, {bt} chips",
        )
    ax.set_xticks(bits)
    ax.set_xlabel("ADC bits (I and Q each)")
    ax.set_ylabel("detection loss against unquantized (dB)")
    ax.set_ylim(0, 6)
    ax.set_xlim(0.7, 4.3)
    ax.legend(loc="upper right", fontsize=label_size)


def hero():
    a2 = read_csv(REPORTS / "adc/a2_losses.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    ax = axes[0]
    _bits_panel(ax)
    ax.set_title("Noise sets the ADC gain: measured loss against theory", fontsize=9.5)

    ax = axes[1]
    series = [
        ("midrise", "none", BLUE, "o", "plain ADC (mid-rise)"),
        ("midrise", "subtractive", ORANGE, "s", "with dither, removed after the ADC"),
        ("midtread", "none", AQUA, "D", "mid-tread ADC (has a zero level)"),
    ]
    for kind, dither, col, mk, label in series:
        pts = sorted(
            (_f(r["a_over_sigma_db"]) - 10 * np.log10(2), _f(r["loss_db"]))
            for r in a2
            if r["bits"] == "3" and r["kind"] == kind and r["dither"] == dither
        )
        x, y = np.array(pts).T
        ax.plot(x, y, color=col, marker=mk, ms=4, lw=1.6, label=label)
    ax.annotate(
        "plain ADC stops at its 1-bit loss",
        (36, 2.34),
        xytext=(0, 10),
        textcoords="offset points",
        ha="center",
        color=INK2,
        fontsize=8,
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.set_ylim(0, 24)
    ax.set_xlim(-16, 44)
    ax.set_xlabel("strong echo elsewhere in the window: SNR per sample (dB)")
    ax.set_ylabel("loss for a weak target (dB)")
    ax.set_title("A strong echo sets the gain (3 bits)", fontsize=9.5)
    fig.suptitle(
        "What a 1 to 4-bit radar receiver costs in detection "
        "(Pd 0.9, Pfa 1e-4, Monte Carlo with 95 % intervals)",
        x=0.01,
        ha="left",
        fontweight="bold",
        fontsize=10.5,
    )
    fig.tight_layout()
    path = ASSETS / "hero.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def social_preview(h):
    fig = plt.figure(figsize=(1280 / 150, 640 / 150), dpi=150, facecolor=SURFACE)
    fig.text(
        0.04,
        0.90,
        "Radar Waveform Design and Matched Filtering:",
        fontsize=15,
        fontweight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.04,
        0.81,
        "LFM, Phase-Coded and OFDM Pulses",
        fontsize=15,
        fontweight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.04, 0.72, "with Low-Bit ADC Analysis", fontsize=14, fontweight="bold", color=INK, va="top"
    )
    fig.text(
        0.04,
        0.60,
        "Detection loss of 1 to 4-bit radar receivers,\n"
        "measured by Monte Carlo and checked against theory",
        fontsize=9.5,
        color=INK2,
        va="top",
    )
    lines = [
        (f"{h['loss_3bit_p4_bt1024_db']['value']:.2f} dB", "3-bit ADC, long code (theory 0.17)"),
        (f"{h['loss_1bit_p4_bt1024_db']['value']:.2f} dB", "1-bit ADC, long code (theory 1.96)"),
        (f"{h['loss_1bit_p4_bt13_db']['value']:.1f} dB", "1-bit ADC, 13-chip code"),
    ]
    y = 0.40
    for big, small in lines:
        fig.text(0.04, y, big, fontsize=16, fontweight="bold", color=BLUE, va="baseline")
        fig.text(0.215, y, small, fontsize=9.5, color=INK2, va="baseline")
        y -= 0.095
    fig.text(0.04, 0.05, "Adan Khalid", fontsize=8.5, color=MUTED)
    ax = fig.add_axes([0.62, 0.17, 0.35, 0.68])
    _bits_panel(ax, label_size=7)
    ax.set_ylabel("loss (dB)")
    ax.set_title("Loss against ADC bits", fontsize=9)
    path = ASSETS / "social-preview.png"
    with plt.rc_context({"savefig.bbox": "standard"}):  # exact 1280 x 640, no trimming
        fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    img = plt.imread(path)
    print(f"wrote {path} ({img.shape[1]} x {img.shape[0]})")


def main():
    h = headline()
    hero()
    social_preview(h)


if __name__ == "__main__":
    main()
