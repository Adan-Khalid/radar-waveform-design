"""Numbers that describe a compressed pulse or a receive filter."""

from __future__ import annotations

import numpy as np

from .matched import correlate


def mainlobe_bounds(mag: np.ndarray, peak: int | None = None) -> tuple[int, int]:
    """Indices of the first minimum on each side of the peak (inclusive bounds)."""
    mag = np.asarray(mag, dtype=float)
    p = int(np.argmax(mag)) if peak is None else int(peak)
    tol = 1e-12 * mag[p]
    left = p
    while left > 0 and mag[left - 1] < mag[left] + tol:
        left -= 1
    right = p
    while right < mag.size - 1 and mag[right + 1] < mag[right] + tol:
        right += 1
    return left, right


def _split(y):
    mag = np.abs(np.asarray(y))
    p = int(np.argmax(mag))
    left, right = mainlobe_bounds(mag, p)
    side = np.concatenate([mag[:left], mag[right + 1 :]])
    return mag, p, mag[left : right + 1], side


def psl_db(y: np.ndarray) -> float:
    """Peak sidelobe level: largest value outside the mainlobe, dB below the peak."""
    mag, p, _, side = _split(y)
    if side.size == 0 or side.max() == 0:
        return float("-inf")
    return float(20 * np.log10(side.max() / mag[p]))


def isl_db(y: np.ndarray) -> float:
    """Integrated sidelobe level: sidelobe energy over mainlobe energy, dB."""
    _, _, main, side = _split(y)
    return float(10 * np.log10(np.sum(side**2) / np.sum(main**2)))


def width(y: np.ndarray, level_db: float = -3.0, samples_per_cell: int = 1) -> float:
    """Mainlobe width at `level_db` below the peak, in range cells (units of 1/B).

    Crossings are found by linear interpolation of |y| between samples.
    """
    mag = np.abs(np.asarray(y))
    p = int(np.argmax(mag))
    level = mag[p] * 10 ** (level_db / 20)

    def crossing(direction):
        i = p
        while 0 <= i + direction < mag.size and mag[i + direction] > level:
            i += direction
        j = i + direction
        if not 0 <= j < mag.size:
            return float(i)
        frac = (mag[i] - level) / (mag[i] - mag[j])
        return i + direction * frac

    return float((crossing(1) - crossing(-1)) / samples_per_cell)


def mismatch_loss_db(s: np.ndarray, ref: np.ndarray) -> float:
    """SNR loss of a receive reference against the matched filter, dB (0 for ref = s).

    loss = -10 log10( max_k |y_k|^2 / (||s||^2 ||ref||^2) ), y = correlate(s, ref).
    By the Cauchy-Schwarz inequality it is never negative.
    """
    y, _ = correlate(s, ref)
    peak = np.max(np.abs(y)) ** 2
    norm = np.sum(np.abs(s) ** 2) * np.sum(np.abs(ref) ** 2)
    return float(-10 * np.log10(peak / norm))


def pulse_metrics(s: np.ndarray, ref: np.ndarray | None = None, samples_per_cell: int = 1):
    """PSL, ISL, -3 dB width (cells) and mismatch loss of s through ref (default: s)."""
    ref = s if ref is None else ref
    y, _ = correlate(s, ref)
    return {
        "psl_db": psl_db(y),
        "isl_db": isl_db(y),
        "width_3db_cells": width(y, -3.0, samples_per_cell),
        "mismatch_loss_db": mismatch_loss_db(s, ref),
    }
