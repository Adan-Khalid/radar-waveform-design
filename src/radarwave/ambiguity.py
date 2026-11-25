"""Narrowband ambiguity function and Doppler tolerance, written from the definition.

chi[m, nu] = sum_n s[n] exp(j 2 pi nu n / f_s) conj(ref[n - m])

It is the output, at lag m, of the filter matched to `ref` when the echo of s arrives
with Doppler shift nu. With ref = s it is the auto-ambiguity function; with a
mismatched receive filter it is the cross-ambiguity function. Doppler is given as nu*T,
cycles of Doppler phase over the pulse, with T = len(s) / f_s.
"""

from __future__ import annotations

import numpy as np

from .metrics import psl_db


def ambiguity(s: np.ndarray, ref: np.ndarray | None = None, nu_t=None, lags=None):
    """Return (chi, lags, nu_t): chi has shape (len(lags), len(nu_t)).

    nu_t=None uses a full DFT grid of K = len(s) points (nu*T = 0, 1, ..., K-1), on which
    the discrete volume identity (1/K) sum |chi|^2 = E_s E_ref holds exactly.
    """
    s = np.asarray(s, dtype=complex)
    ref = s if ref is None else np.asarray(ref, dtype=complex)
    n_s, n_r = s.size, ref.size
    lags = np.arange(-(n_r - 1), n_s) if lags is None else np.asarray(lags)
    nu_t = np.arange(n_s, dtype=float) if nu_t is None else np.asarray(nu_t, dtype=float)
    n = np.arange(n_s)
    idx = n[None, :] - lags[:, None]
    valid = (idx >= 0) & (idx < n_r)
    prod = np.where(valid, s[None, :] * np.conj(ref[np.clip(idx, 0, n_r - 1)]), 0)
    doppler = np.exp(2j * np.pi * np.outer(n, nu_t) / n_s)
    return prod @ doppler, lags, nu_t


def doppler_cut_summary(s, ref=None, nu_t=None, samples_per_cell: int = 1, chunk: int = 256):
    """Peak loss, delay shift and PSL of the compressed pulse at each Doppler value.

    Returns a dict of arrays over nu_t:
    - peak_loss_db: 20 log10(max over lag |chi(., nu)| / max over lag |chi(., 0)|)
    - shift_cells: lag of the peak, in range cells
    - psl_db: peak sidelobe relative to the peak at that Doppler
    """
    nu_t = np.asarray(nu_t, dtype=float)
    peak = np.empty(nu_t.size)
    shift = np.empty(nu_t.size)
    psl = np.empty(nu_t.size)
    for start in range(0, nu_t.size, chunk):
        sl = slice(start, start + chunk)
        chi, lags, _ = ambiguity(s, ref, nu_t[sl])
        mag = np.abs(chi)
        k = np.argmax(mag, axis=0)
        peak[sl] = mag[k, np.arange(mag.shape[1])]
        shift[sl] = lags[k] / samples_per_cell
        psl[sl] = [psl_db(mag[:, j]) for j in range(mag.shape[1])]
    chi0, _, _ = ambiguity(s, ref, [0.0])
    peak0 = np.max(np.abs(chi0))
    return {
        "nu_t": nu_t,
        "peak_loss_db": 20 * np.log10(peak / peak0),
        "shift_cells": shift,
        "psl_db": psl,
    }


def first_crossing(x: np.ndarray, y: np.ndarray, level: float) -> float:
    """Smallest x (x ascending, from 0) where y first falls below `level`; linear interpolation."""
    below = np.nonzero(y < level)[0]
    if below.size == 0:
        return float("inf")
    j = below[0]
    if j == 0:
        return float(x[0])
    x0, x1, y0, y1 = x[j - 1], x[j], y[j - 1], y[j]
    return float(x0 + (level - y0) * (x1 - x0) / (y1 - y0))


def doppler_tolerance(s, ref=None, nu_t_max: float = 4.0, step: float = 0.01, samples_per_cell=1):
    """Summary numbers for one waveform, from nu*T = 0 to nu_t_max:

    - nu_t_1db, nu_t_3db: first nu*T where the peak has lost 1 dB and 3 dB (inf if never)
    - loss_at_<v>_db for nu*T = 0.5, 1, 2, 4: peak loss there. Polyphase codes (P4,
      Frank) dip at half-integer nu*T and recover near integer nu*T, one chip further on,
      so the first crossing alone undersells them.
    - shift_at_<v>_cells for nu*T = 1, 4: delay of the peak in range cells
      (range-Doppler coupling: LFM moves about one cell per unit of nu*T)
    - psl_at_0p25_db, psl_at_0p5_db: PSL at nu*T = 0.25 and 0.5
    """
    nu = np.round(np.arange(0.0, nu_t_max + step / 2, step), 10)
    cut = doppler_cut_summary(s, ref, nu, samples_per_cell)

    def at(key, v):
        j = int(np.argmin(np.abs(nu - v)))
        return float(cut[key][j]) if abs(nu[j] - v) < step / 2 else float("nan")

    out = {
        "nu_t_1db": first_crossing(nu, cut["peak_loss_db"], -1.0),
        "nu_t_3db": first_crossing(nu, cut["peak_loss_db"], -3.0),
    }
    for v, tag in ((0.5, "0p5"), (1.0, "1"), (2.0, "2"), (4.0, "4")):
        out[f"loss_at_{tag}_db"] = at("peak_loss_db", v)
    for v, tag in ((1.0, "1"), (4.0, "4")):
        out[f"shift_at_{tag}_cells"] = at("shift_cells", v)
    out["psl_at_0p25_db"] = at("psl_db", 0.25)
    out["psl_at_0p5_db"] = at("psl_db", 0.5)
    return out, cut
