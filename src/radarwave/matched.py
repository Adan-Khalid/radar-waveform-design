"""Matched and mismatched filtering by FFT cross-correlation.

Convention: correlate(x, ref)[k] = sum_n x[n + k] conj(ref[n]). An echo of `ref` that
starts at sample d of x gives the matched peak at lag k = d.
"""

from __future__ import annotations

import numpy as np


def _nfft(n: int) -> int:
    return 1 << int(np.ceil(np.log2(max(n, 1))))


def correlate(x: np.ndarray, ref: np.ndarray, lags: str = "full"):
    """Cross-correlate x (last axis) with a reference.

    lags="full" returns (y, lag) for lags -(len(ref)-1) .. len(x)-1.
    lags="valid" returns (y, lag) for lags 0 .. len(x)-len(ref), where the reference
    lies entirely inside x (every output uses the same number of samples).
    """
    x = np.asarray(x)
    ref = np.asarray(ref)
    m, n_ref = x.shape[-1], ref.size
    nfft = _nfft(m + n_ref - 1)
    spec = np.fft.fft(x, nfft, axis=-1) * np.conj(np.fft.fft(ref, nfft))
    c = np.fft.ifft(spec, axis=-1)
    if lags == "valid":
        if m < n_ref:
            raise ValueError("x is shorter than the reference")
        return c[..., : m - n_ref + 1], np.arange(m - n_ref + 1)
    if lags != "full":
        raise ValueError("lags must be 'full' or 'valid'")
    y = np.concatenate([c[..., nfft - n_ref + 1 :], c[..., :m]], axis=-1)
    return y, np.arange(-(n_ref - 1), m)


def matched_filter(x: np.ndarray, s: np.ndarray):
    """Matched filter output of x for the waveform s (the reference is s itself)."""
    return correlate(x, s)


def autocorrelation(s: np.ndarray):
    """Compressed pulse of s through its own matched filter: lags -(N-1)..N-1."""
    return correlate(s, s)
