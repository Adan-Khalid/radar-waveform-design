"""Closed-form and numerically integrated predictions that the simulations are checked against."""

from __future__ import annotations

import numpy as np
from scipy import optimize, special, stats

from .adc import staircase, step


def pd_marcum(snr_db, pfa: float):
    """Pd of a non-fluctuating target with random phase, square-law detector, one pulse.

    Under noise only |y|^2 / sigma^2 is exponential, so the threshold is -ln(Pfa).
    With signal, 2|y|^2/sigma^2 is non-central chi-square (2 degrees of freedom,
    non-centrality 2 SNR): Pd = Q1(sqrt(2 SNR), sqrt(-2 ln Pfa)), Marcum's Q function.
    """
    snr = 10 ** (np.asarray(snr_db, dtype=float) / 10)
    return stats.ncx2.sf(-2 * np.log(pfa), 2, 2 * snr)


def required_snr_db(pd: float, pfa: float) -> float:
    """SNR (dB) at which pd_marcum reaches pd."""
    return float(optimize.brentq(lambda x: pd_marcum(x, pfa) - pd, -20, 40, xtol=1e-10))


def albersheim_db(pd: float, pfa: float, n_pulses: int = 1) -> float:
    """Albersheim's approximation (linear detector, non-fluctuating), for comparison."""
    a = np.log(0.62 / pfa)
    b = np.log(pd / (1 - pd))
    return float(
        -5 * np.log10(n_pulses)
        + (6.2 + 4.54 / np.sqrt(n_pulses + 0.44)) * np.log10(a + 0.12 * a * b + 1.7 * b)
    )


def window_loss_db(w: np.ndarray) -> float:
    """SNR loss of an amplitude weighting on a constant-envelope pulse: (mean w)^2/mean(w^2)."""
    w = np.asarray(w, dtype=float)
    return float(-10 * np.log10(np.mean(w) ** 2 / np.mean(w**2)))


def _moments(thresholds, levels, shift):
    """E[x q(x+shift)], E[q(x+shift)], E[q(x+shift)^2] for x ~ N(0, 1), exactly.

    `shift` may be an array; the results then have its shape.
    """
    shift = np.asarray(shift, dtype=float)[..., None]
    lo = np.concatenate([[-np.inf], thresholds]) - shift
    hi = np.concatenate([thresholds, [np.inf]]) - shift
    prob = special.ndtr(hi) - special.ndtr(lo)
    pdf_lo = np.exp(-0.5 * lo**2) / np.sqrt(2 * np.pi)  # exp(-inf) = 0 at the open ends
    pdf_hi = np.exp(-0.5 * hi**2) / np.sqrt(2 * np.pi)
    exq = np.sum(levels * (pdf_lo - pdf_hi), axis=-1)
    eq = np.sum(levels * prob, axis=-1)
    eq2 = np.sum(levels**2 * prob, axis=-1)
    return exq, eq, eq2


def quantizer_efficiency(
    bits: int | None,
    full_scale_over_sigma: float,
    kind: str = "midrise",
    dither: str = "none",
    nodes: int = 64,
) -> float:
    """Low-SNR efficiency of a quantizer on a Gaussian input (Bussgang).

        eta = E[x z]^2 / (sigma^2 E[z^2]),

    with z the receiver's output for input x (after dither and, for subtractive dither,
    after removing it). A weak signal in Gaussian noise keeps a fraction eta of its SNR,
    so the predicted loss is -10 log10(eta). Expectations over x are exact (normal CDF
    sums); the average over the uniform dither uses Gauss-Legendre nodes.
    full_scale_over_sigma: A / sigma of one component.
    """
    if bits is None:
        return 1.0
    thr, lev = staircase(bits, full_scale_over_sigma, kind)
    if dither == "none":
        exq, _, eq2 = _moments(thr, lev, 0.0)
        if eq2 <= 0.0:  # a mid-tread quantizer far too coarse outputs only zeros
            return 0.0
        return float(exq**2 / eq2)
    d = step(bits, full_scale_over_sigma)
    g, wts = np.polynomial.legendre.leggauss(nodes)
    shifts = g * d / 2  # dither values in [-d/2, d/2]
    wts = wts / 2  # average over the interval
    exq, eq, eq2 = _moments(thr, lev, shifts)
    exz = np.sum(wts * exq)
    if dither == "nonsubtractive":
        ez2 = np.sum(wts * eq2)
    elif dither == "subtractive":
        ez2 = np.sum(wts * (eq2 - 2 * shifts * eq + shifts**2))
    else:
        raise ValueError("dither must be none, nonsubtractive or subtractive")
    return float(exz**2 / ez2)


def quantizer_loss_db(bits, full_scale_over_sigma, kind="midrise", dither="none") -> float:
    """-10 log10 of quantizer_efficiency: predicted SNR loss in dB at low per-sample SNR."""
    eff = quantizer_efficiency(bits, full_scale_over_sigma, kind, dither)
    return float("inf") if eff <= 0 else float(-10 * np.log10(eff))


def lfm_ridge_peak(nu_t: np.ndarray, time_bandwidth: float) -> np.ndarray:
    """Height of the LFM ambiguity ridge at Doppler nu: 1 - |nu| / B = 1 - |nu T| / BT."""
    return np.clip(1 - np.abs(np.asarray(nu_t)) / time_bandwidth, 0, None)


def dirichlet_zero_delay(nu_t: np.ndarray, n_samples: int) -> np.ndarray:
    """|chi(0, nu)| of a constant-envelope pulse of n samples (a sampled sinc)."""
    nu_t = np.atleast_1d(np.asarray(nu_t, dtype=float))
    num = np.sin(np.pi * nu_t)
    den = n_samples * np.sin(np.pi * nu_t / n_samples)
    out = np.ones_like(nu_t)
    nz = np.abs(den) > 1e-15
    out[nz] = np.abs(num[nz] / den[nz])
    return out
