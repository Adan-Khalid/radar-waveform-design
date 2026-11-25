"""Monte Carlo detection through a quantizing receiver, with confidence intervals.

Model (one complex sample per range cell, f_s = B):
    r[n] = a exp(j phi) s[n] + w[n],   w ~ CN(0, 1),   sum |s|^2 = 1,
so the unquantized matched-filter SNR is |a|^2 (`snr_db`). Each receiver quantizes I and
Q (`adc.quantize_with_dither`) with its full scale given in units of the noise standard
deviation of one component, sigma_c = 1/sqrt(2). The statistic is |sum_n z[n] s*[n]|^2 at
the target's cell; the threshold of each receiver comes from its own noise-only
statistics at the chosen Pfa (clairvoyant detector: no CFAR loss, no straddle loss).

Common random numbers: every receiver sees the same noise, target phases and dither
draws, and the same noise is reused at every SNR (only the echo amplitude changes), so
differences between receivers are measured with little Monte Carlo noise and can be
bootstrapped in pairs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import theory
from .adc import gaussian_optimal_full_scale, quantize_with_dither
from .matched import correlate

SIGMA_C = 1 / np.sqrt(2)  # noise standard deviation of I or Q (complex noise power 1)


@dataclass(frozen=True)
class Receiver:
    """ADC settings of one receiver. full_scale: A / sigma_c (None: Max's Gaussian optimum)."""

    bits: int | None = None
    kind: str = "midrise"
    dither: str = "none"
    full_scale: float | None = None

    @property
    def a_over_sigma(self) -> float:
        if self.bits is None:
            return float("nan")
        if self.full_scale is None:
            return gaussian_optimal_full_scale(self.bits)
        return float(self.full_scale)

    @property
    def label(self) -> str:
        if self.bits is None:
            return "unquantized"
        parts = [f"{self.bits} bit"]
        if self.kind != "midrise":
            parts.append(self.kind)
        if self.dither != "none":
            parts.append(f"{self.dither} dither")
        if self.full_scale is not None:
            parts.append(f"A/sigma {self.full_scale:.3g}")
        return ", ".join(parts)

    def theory_loss_db(self) -> float:
        """Low-SNR prediction (Bussgang efficiency) of this receiver's SNR loss."""
        if self.bits is None:
            return 0.0
        return theory.quantizer_loss_db(self.bits, self.a_over_sigma, self.kind, self.dither)

    def apply(self, re: np.ndarray, im: np.ndarray, u_re=None, u_im=None) -> np.ndarray:
        """Quantize I and Q (float arrays) and return complex samples."""
        if self.bits is None:
            return re + 1j * im
        a = self.a_over_sigma * SIGMA_C
        qi = quantize_with_dither(re, self.bits, a, self.kind, self.dither, u_re)
        qq = quantize_with_dither(im, self.bits, a, self.kind, self.dither, u_im)
        return qi + 1j * qq


# ---------------------------------------------------------------- intervals


def wilson(k, n, z: float = 1.959964):
    """Wilson score interval for a binomial proportion k/n (arrays allowed)."""
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    p = k / n
    den = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return centre - half, centre + half


def crossing_db(snr_db: np.ndarray, pd: np.ndarray, target: float) -> np.ndarray:
    """SNR where Pd reaches `target`, by linear interpolation above the last point below it.

    pd may have a trailing axis of bootstrap replicates: shape (n_snr,) or (n_snr, R).
    Returns inf where Pd never reaches the target and -inf where it starts above it.
    """
    pd = np.asarray(pd, dtype=float)
    squeeze = pd.ndim == 1
    if squeeze:
        pd = pd[:, None]
    n_snr, n_rep = pd.shape
    out = np.empty(n_rep)
    below = pd < target
    for r in range(n_rep):
        idx = np.nonzero(below[:, r])[0]
        if idx.size == 0:
            out[r] = -np.inf
            continue
        j = idx[-1]
        if j == n_snr - 1:
            out[r] = np.inf
            continue
        p0, p1 = pd[j, r], pd[j + 1, r]
        out[r] = snr_db[j] + (target - p0) * (snr_db[j + 1] - snr_db[j]) / (p1 - p0)
    return out[0] if squeeze else out


# ---------------------------------------------------------------- noise-only statistics


def noise_statistics(
    s: np.ndarray,
    receivers: list[Receiver],
    n_stats: int,
    seed,
    record_factor: int = 16,
    rows_per_batch: int | None = None,
) -> dict[Receiver, np.ndarray]:
    """Noise-only detection statistics of every receiver (same noise for all).

    Long noise records are quantized once and correlated at every lag where the
    waveform fits inside the record, so each record gives many statistics.
    """
    s = np.asarray(s, dtype=np.complex64)
    n = s.size
    m = record_factor * n
    per_row = m - n + 1
    rows_total = int(np.ceil(n_stats / per_row))
    rows_per_batch = rows_per_batch or max(1, min(rows_total, 4_000_000 // m))
    out = {rx: [] for rx in receivers}
    done = 0
    batch = 0
    while done < rows_total:
        rows = min(rows_per_batch, rows_total - done)
        rng = np.random.default_rng([*np.atleast_1d(seed), batch])
        re = rng.standard_normal((rows, m), dtype=np.float32) * np.float32(SIGMA_C)
        im = rng.standard_normal((rows, m), dtype=np.float32) * np.float32(SIGMA_C)
        u_re = rng.random((rows, m), dtype=np.float32) - np.float32(0.5)
        u_im = rng.random((rows, m), dtype=np.float32) - np.float32(0.5)
        for rx in receivers:
            z = rx.apply(re, im, u_re, u_im)
            y, _ = correlate(z, s, "valid")
            out[rx].append((np.abs(y) ** 2).astype(np.float32).ravel())
        done += rows
        batch += 1
    return {rx: np.concatenate(v)[:n_stats] for rx, v in out.items()}


TIE_TOL = 1e-4  # relative width inside which two statistics count as the same value


@dataclass(frozen=True)
class Threshold:
    """Neyman-Pearson threshold, randomised on ties.

    After a coarse quantizer the statistic can take only a few values (one bit and a
    short binary code give a handful), so no plain threshold gives Pfa exactly. The
    test decides "target" when stat > value, and with probability gamma when stat
    equals value; gamma is set from the noise-only statistics so that the false-alarm
    rate is the design value. For continuous statistics gamma has no effect.
    """

    value: float
    gamma: float = 0.0

    def decide(self, stat: np.ndarray, u: np.ndarray | None = None) -> np.ndarray:
        band = TIE_TOL * self.value
        above = stat > self.value + band
        if self.gamma <= 0.0 or u is None:
            return above
        tie = np.abs(stat - self.value) <= band
        return above | (tie & (u < self.gamma))


def thresholds_for(stats: dict[Receiver, np.ndarray], pfa: float) -> dict[Receiver, Threshold]:
    """Randomised thresholds at the design Pfa from noise-only statistics."""
    out = {}
    for rx, v in stats.items():
        t = float(np.quantile(v, 1 - pfa, method="higher"))  # an observed value
        band = TIE_TOL * t
        n = v.size
        n_above = int(np.sum(v > t + band))
        n_tie = int(np.sum(np.abs(v - t) <= band))
        gamma = float(np.clip((pfa * n - n_above) / n_tie, 0.0, 1.0)) if n_tie else 0.0
        out[rx] = Threshold(t, gamma)
    return out


def check_pfa(s, receivers, thresholds, n_stats, seed, pfa):
    """False-alarm rate of each threshold on fresh noise, with Wilson intervals."""
    stats = noise_statistics(s, receivers, n_stats, seed)
    rng = np.random.default_rng([*np.atleast_1d(seed), 31337])
    rows = []
    for rx in receivers:
        u = rng.random(stats[rx].size)
        k = int(np.sum(thresholds[rx].decide(stats[rx], u)))
        n = stats[rx].size
        lo, hi = wilson(k, n)
        rows.append(
            {
                "receiver": rx.label,
                "false_alarms": k,
                "statistics": n,
                "pfa_design": pfa,
                "pfa_measured": k / n,
                "pfa_lo": float(lo),
                "pfa_hi": float(hi),
                "design_inside": bool(lo <= pfa <= hi),
            }
        )
    return rows


# ---------------------------------------------------------------- detections


class TrialSet:
    """Fixed noise, phases and dither draws for detection trials (common random numbers)."""

    def __init__(self, n_samples: int, trials: int, seed, batch: int = 4000):
        self.n = n_samples
        self.trials = trials
        self.seed = seed
        self.batch = batch

    def batches(self):
        for b, start in enumerate(range(0, self.trials, self.batch)):
            rows = min(self.batch, self.trials - start)
            rng = np.random.default_rng([*np.atleast_1d(self.seed), 7919, b])
            sc = np.float32(SIGMA_C)
            re = rng.standard_normal((rows, self.n), dtype=np.float32) * sc
            im = rng.standard_normal((rows, self.n), dtype=np.float32) * sc
            phase = rng.random(rows) * 2 * np.pi
            u_re = rng.random((rows, self.n), dtype=np.float32) - np.float32(0.5)
            u_im = rng.random((rows, self.n), dtype=np.float32) - np.float32(0.5)
            u_tie = rng.random(rows)
            yield slice(start, start + rows), re, im, phase, u_re, u_im, u_tie


def detections(
    s, rx: Receiver, threshold: Threshold, snr_db, trials: TrialSet, phase_deg: float | None = None
) -> np.ndarray:
    """Boolean detections, shape (len(snr_db), trials), at the target's cell.

    phase_deg fixes the echo phase for every trial (default: random phase per trial).
    """
    s = np.asarray(s, dtype=np.complex64)
    snr_db = np.atleast_1d(np.asarray(snr_db, dtype=float))
    hits = np.zeros((snr_db.size, trials.trials), dtype=bool)
    s_conj = np.conj(s)
    if not isinstance(threshold, Threshold):
        threshold = Threshold(float(threshold))
    for sl, re, im, phase, u_re, u_im, u_tie in trials.batches():
        if phase_deg is not None:
            phase = np.full(phase.shape, np.deg2rad(phase_deg))
        echo = np.exp(1j * phase)[:, None] * s[None, :]
        e_re = echo.real.astype(np.float32)
        e_im = echo.imag.astype(np.float32)
        for j, snr in enumerate(snr_db):
            amp = np.float32(10 ** (snr / 20))
            x_re = re + amp * e_re
            x_im = im + amp * e_im
            z = rx.apply(x_re, x_im, u_re, u_im)
            stat = np.abs(z @ s_conj) ** 2
            hits[j, sl] = threshold.decide(stat, u_tie)
    return hits


def poisson_weights(trials: int, replicates: int, seed) -> np.ndarray:
    """Poisson(1) bootstrap weights, shape (trials, replicates), float32."""
    rng = np.random.default_rng([*np.atleast_1d(seed), 104729])
    return rng.poisson(1.0, size=(trials, replicates)).astype(np.float32)


@dataclass
class PdCurve:
    receiver: Receiver
    snr_db: np.ndarray
    hits: np.ndarray  # (n_snr, trials) bool

    @property
    def pd(self) -> np.ndarray:
        return self.hits.mean(axis=1)

    def required_snr(self, target: float) -> float:
        return float(crossing_db(self.snr_db, self.pd, target))

    def required_snr_bootstrap(self, target: float, weights: np.ndarray) -> np.ndarray:
        pd_rep = (self.hits.astype(np.float32) @ weights) / weights.sum(axis=0)
        return crossing_db(self.snr_db, pd_rep, target)


def search_required_snr(
    s,
    rx: Receiver,
    threshold: Threshold,
    trials: TrialSet,
    coarse_trials: TrialSet,
    target: float,
    predicted_db: float,
    fine_half_width: float = 1.0,
    fine_step: float = 0.1,
    coarse_step: float = 0.5,
    max_snr_db: float = 60.0,
) -> PdCurve:
    """Coarse grid with few trials to bracket Pd = target, then a fine grid with all trials."""
    lo = predicted_db - 3.0
    hi = min(predicted_db + 9.0, max_snr_db)
    approx = None
    while approx is None:
        grid = np.arange(lo, hi + coarse_step / 2, coarse_step)
        pd = detections(s, rx, threshold, grid, coarse_trials).mean(axis=1)
        x = crossing_db(grid, pd, target)
        if np.isfinite(x):
            approx = x
        elif x == -np.inf:  # already above the target at the bottom of the grid
            hi, lo = lo + coarse_step, lo - 6.0
        elif hi >= max_snr_db:  # never reached
            approx = np.inf
        else:
            lo, hi = hi - coarse_step, min(hi + 9.0, max_snr_db)
    if not np.isfinite(approx):
        grid = np.arange(max_snr_db - 2 * fine_half_width, max_snr_db + fine_step / 2, fine_step)
        return PdCurve(rx, grid, detections(s, rx, threshold, grid, trials))
    grid = np.round(
        np.arange(approx - fine_half_width, approx + fine_half_width + fine_step / 2, fine_step),
        6,
    )
    hits = detections(s, rx, threshold, grid, trials)
    curve = PdCurve(rx, grid, hits)
    x = curve.required_snr(target)
    if not np.isfinite(x):  # widen once if the fine grid missed the crossing
        grid2 = np.round(np.arange(approx - 3, approx + 3 + fine_step / 2, fine_step), 6)
        curve = PdCurve(rx, grid2, detections(s, rx, threshold, grid2, trials))
    return curve


def loss_table(
    curves: dict[Receiver, PdCurve],
    reference: Receiver,
    target: float,
    weights: np.ndarray,
    pfa: float,
) -> list[dict]:
    """Required SNR and loss against the reference, with paired bootstrap intervals."""
    ref_boot = curves[reference].required_snr_bootstrap(target, weights)
    ref_req = curves[reference].required_snr(target)
    req_theory = theory.required_snr_db(target, pfa)
    rows = []
    for rx, curve in curves.items():
        req = curve.required_snr(target)
        boot = curve.required_snr_bootstrap(target, weights)
        diff = boot - ref_boot
        finite = np.isfinite(boot)
        req_lo, req_hi = np.percentile(boot[finite], [2.5, 97.5]) if finite.all() else (np.nan,) * 2
        loss_lo, loss_hi = (
            np.percentile(diff[finite], [2.5, 97.5]) if finite.all() else (np.nan,) * 2
        )
        rows.append(
            {
                "receiver": rx.label,
                "bits": "" if rx.bits is None else rx.bits,
                "kind": rx.kind,
                "dither": rx.dither,
                "a_over_sigma": rx.a_over_sigma,
                "required_snr_db": req,
                "required_snr_lo": float(req_lo),
                "required_snr_hi": float(req_hi),
                "loss_db": req - ref_req,
                "loss_lo": float(loss_lo),
                "loss_hi": float(loss_hi),
                "theory_loss_db": rx.theory_loss_db(),
                "theory_required_snr_db": req_theory + rx.theory_loss_db(),
            }
        )
    return rows


def run_receivers(
    s,
    receivers: list[Receiver],
    *,
    pfa: float,
    pd_target: float,
    noise_stats: int,
    trials: int,
    coarse_trials: int,
    bootstrap: int,
    seed: int,
    check_stats: int | None = None,
    fine_half_width: float = 1.0,
    fine_step: float = 0.1,
):
    """Thresholds, Pfa check, Pd curves and paired losses for a set of receivers.

    The first receiver in the list is the reference for the losses (normally unquantized).
    """
    stats = noise_statistics(s, receivers, noise_stats, [seed, 1])
    thr = thresholds_for(stats, pfa)
    pfa_rows = check_pfa(s, receivers, thr, check_stats or noise_stats // 2, [seed, 2], pfa)
    trial_set = TrialSet(np.size(s), trials, [seed, 3])
    coarse_set = TrialSet(np.size(s), coarse_trials, [seed, 4])
    req0 = theory.required_snr_db(pd_target, pfa)
    curves = {}
    for rx in receivers:
        pred = req0 + min(rx.theory_loss_db(), 40.0)
        curves[rx] = search_required_snr(
            s,
            rx,
            thr[rx],
            trial_set,
            coarse_set,
            pd_target,
            pred,
            fine_half_width=fine_half_width,
            fine_step=fine_step,
        )
    weights = poisson_weights(trials, bootstrap, [seed, 5])
    loss_rows = loss_table(curves, receivers[0], pd_target, weights, pfa)
    for row, rx in zip(loss_rows, receivers, strict=True):
        row["threshold"] = thr[rx].value
        row["threshold_gamma"] = thr[rx].gamma
    curve_rows = []
    for rx, c in curves.items():
        lo, hi = wilson(c.hits.sum(axis=1), c.hits.shape[1])
        for j, snr in enumerate(c.snr_db):
            curve_rows.append(
                {
                    "receiver": rx.label,
                    "snr_db": float(snr),
                    "pd": float(c.pd[j]),
                    "pd_lo": float(lo[j]),
                    "pd_hi": float(hi[j]),
                    "trials": c.hits.shape[1],
                    "pd_marcum": float(theory.pd_marcum(snr, pfa)),
                }
            )
    return loss_rows, curve_rows, pfa_rows, curves
