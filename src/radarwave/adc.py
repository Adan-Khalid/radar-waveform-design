"""ADC model: gain (full scale), uniform quantizer of 1 to 4 bits, optional dither.

Each of I and Q is quantized on its own. The full scale A is the clipping level of one
component; the step is 2A / 2^bits for both quantizer kinds:

- mid-rise (default): 2^bits levels at +-step/2, +-3 step/2, ..., no zero level.
  One bit keeps only the sign.
- mid-tread: 2^bits - 1 levels at 0, +-step, ..., +-(2^(bits-1) - 1) step (needs 2+ bits).

Dither is uniform over one step (from -step/2 to +step/2), independent per sample and
per component. Non-subtractive: quantize(x + d). Subtractive: quantize(x + d) - d, as a
receiver with a known pseudo-random dither source would do.
"""

from __future__ import annotations

import numpy as np

# Max (1960), Table II: best uniform step for a zero-mean, unit-variance Gaussian input,
# for 2, 4, 8 and 16 levels (1, 2, 3 and 4 bits).
GAUSSIAN_OPTIMAL_STEP = {1: 1.596, 2: 0.9957, 3: 0.5860, 4: 0.3352}

KINDS = ("midrise", "midtread")
DITHERS = ("none", "nonsubtractive", "subtractive")


def step(bits: int, full_scale: float) -> float:
    """Quantizer step for a full scale A (clipping level of one component)."""
    return 2.0 * full_scale / 2**bits


def gaussian_optimal_full_scale(bits: int, sigma: float = 1.0) -> float:
    """Full scale A that gives Max's optimal mid-rise step for a Gaussian of std sigma."""
    return 2 ** (bits - 1) * GAUSSIAN_OPTIMAL_STEP[bits] * sigma


def quantize_real(x: np.ndarray, bits: int | None, full_scale: float, kind: str = "midrise"):
    """Uniform quantizer on real samples; bits=None returns x unchanged."""
    if bits is None:
        return x
    d = step(bits, full_scale)
    if kind == "midrise":
        half = 2 ** (bits - 1)
        k = np.clip(np.floor(x / d), -half, half - 1)
        return (k + 0.5) * d
    if kind == "midtread":
        if bits < 2:
            raise ValueError("a mid-tread quantizer needs at least 2 bits")
        kmax = 2 ** (bits - 1) - 1
        k = np.clip(np.floor(x / d + 0.5), -kmax, kmax)
        return k * d
    raise ValueError(f"kind must be one of {KINDS}")


def staircase(bits: int, full_scale: float, kind: str = "midrise"):
    """Decision thresholds (ascending) and output levels of the quantizer."""
    d = step(bits, full_scale)
    if kind == "midrise":
        half = 2 ** (bits - 1)
        levels = (np.arange(-half, half) + 0.5) * d
        thresholds = np.arange(-half + 1, half) * d
    elif kind == "midtread":
        if bits < 2:
            raise ValueError("a mid-tread quantizer needs at least 2 bits")
        kmax = 2 ** (bits - 1) - 1
        levels = np.arange(-kmax, kmax + 1) * d
        thresholds = (np.arange(-kmax, kmax) + 0.5) * d
    else:
        raise ValueError(f"kind must be one of {KINDS}")
    return thresholds, levels


def quantize_with_dither(
    x: np.ndarray,
    bits: int | None,
    full_scale: float,
    kind: str = "midrise",
    dither: str = "none",
    uniform: np.ndarray | None = None,
) -> np.ndarray:
    """Quantize real samples, with dither drawn as `uniform` (values in [-1/2, 1/2))."""
    if bits is None or dither == "none":
        return quantize_real(x, bits, full_scale, kind)
    if uniform is None:
        raise ValueError("dither needs an array of uniform values in [-1/2, 1/2)")
    dv = uniform * step(bits, full_scale)
    q = quantize_real(x + dv, bits, full_scale, kind)
    if dither == "nonsubtractive":
        return q
    if dither == "subtractive":
        return q - dv
    raise ValueError(f"dither must be one of {DITHERS}")


def quantize(
    z: np.ndarray,
    bits: int | None,
    full_scale: float,
    kind: str = "midrise",
    dither: str = "none",
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Quantize complex samples: I and Q each through the same quantizer."""
    z = np.asarray(z)
    if bits is None:
        return z
    if dither == "none":
        u_i = u_q = None
    else:
        rng = rng if rng is not None else np.random.default_rng()
        u_i = rng.random(z.shape) - 0.5
        u_q = rng.random(z.shape) - 0.5
    i = quantize_with_dither(z.real, bits, full_scale, kind, dither, u_i)
    q = quantize_with_dither(z.imag, bits, full_scale, kind, dither, u_q)
    return i + 1j * q
