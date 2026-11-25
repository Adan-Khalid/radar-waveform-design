"""Sidelobe control: amplitude windows on the receive reference, and mismatched filters.

Every function here returns a receive *reference* r for `matched.correlate(x, r)`, the
same form as the matched filter (whose reference is the waveform itself).
"""

from __future__ import annotations

import numpy as np
from scipy import optimize


def generalized_cosine(n: int, a: tuple[float, ...]) -> np.ndarray:
    """Symmetric window sum_k (-1)^k a_k cos(2 pi k i / (n - 1)), i = 0..n-1."""
    if n == 1:
        return np.ones(1)
    i = np.arange(n)
    return sum((-1) ** k * ak * np.cos(2 * np.pi * k * i / (n - 1)) for k, ak in enumerate(a))


def taylor(n: int, sidelobe_db: float = 35.0, nbar: int = 5) -> np.ndarray:
    """Taylor window (the formula MATLAB's taylorwin uses; equals SciPy's with norm=False)."""
    a = np.arccosh(10 ** (sidelobe_db / 20)) / np.pi
    sp2 = nbar**2 / (a**2 + (nbar - 0.5) ** 2)
    m = np.arange(1, nbar)
    f = np.empty(nbar - 1)
    for j, mj in enumerate(m):
        numer = np.prod(1 - mj**2 / sp2 / (a**2 + (m - 0.5) ** 2))
        others = m[m != mj]
        denom = 2 * np.prod(1 - mj**2 / others**2)
        f[j] = (-1) ** (mj + 1) * numer / denom
    x = (np.arange(n) - n / 2 + 0.5) / n
    return 1 + 2 * np.cos(2 * np.pi * np.outer(x, m)) @ f


def window(name: str, n: int, **kw) -> np.ndarray:
    """Symmetric amplitude window by name."""
    if name == "rectangular":
        return np.ones(n)
    if name == "hann":
        return generalized_cosine(n, (0.5, 0.5))
    if name == "hamming":
        return generalized_cosine(n, (0.54, 0.46))
    if name == "blackman":
        return generalized_cosine(n, (0.42, 0.5, 0.08))
    if name == "taylor":
        return taylor(n, kw.get("sidelobe_db", 35.0), kw.get("nbar", 5))
    raise ValueError(f"unknown window {name!r}")


def window_label(spec: dict) -> str:
    name = spec["name"]
    if name == "taylor":
        return f"Taylor {spec.get('sidelobe_db', 35)} dB, nbar {spec.get('nbar', 5)}"
    return name.capitalize()


def windowed_reference(s: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Receive reference with amplitude weighting w applied over the pulse."""
    return np.asarray(s) * np.asarray(w)


def convolution_matrix(s: np.ndarray, length: int) -> np.ndarray:
    """Matrix S with S @ g = full convolution of s with g (g of `length` taps)."""
    s = np.asarray(s)
    n_out = s.size + length - 1
    i = np.arange(n_out)[:, None] - np.arange(length)[None, :]
    valid = (i >= 0) & (i < s.size)
    return np.where(valid, s[np.clip(i, 0, s.size - 1)], 0)


def _to_reference(g: np.ndarray) -> np.ndarray:
    # correlate(x, r) with r = conj(reversed g) equals convolution of x with g.
    return np.conj(g[::-1])


def ls_mismatched_filter(s: np.ndarray, length: int, loading: float = 1e-9) -> np.ndarray:
    """Least-squares mismatched filter (minimum ISL with the peak fixed), as a reference.

    With z = S g the output (S: convolution matrix of the code, g: filter taps), fix
    z[c] = 1 at the centre c and minimise the energy of all other outputs:
        g = R^-1 a^H / (a R^-1 a^H),  R = S_rest^H S_rest + loading * I,
    a = row c of S, S_rest = S without row c (Ackroyd and Ghani's form).
    """
    s = np.asarray(s, dtype=complex)
    big = convolution_matrix(s, length)
    c = (big.shape[0] - 1) // 2
    a = big[c]
    rest = np.delete(big, c, axis=0)
    r = rest.conj().T @ rest + loading * np.sum(np.abs(s) ** 2) * np.eye(length)
    g = np.linalg.solve(r, a.conj())
    g = g / (a @ g)
    return _to_reference(g)


def minpsl_filter_real(s: np.ndarray, length: int) -> np.ndarray:
    """Minimum peak sidelobe filter for a real (binary) code, by linear programming.

    minimise t subject to |z_k| <= t for k != c and z_c = 1, z = S g, g real
    (Baden and Cohen's criterion). Returned as a reference.
    """
    s = np.asarray(s)
    if np.iscomplexobj(s) and np.any(np.abs(s.imag) > 1e-12):
        raise ValueError("the linear-programming filter here is for real codes")
    big = convolution_matrix(s.real, length).real
    c = (big.shape[0] - 1) // 2
    rest = np.delete(big, c, axis=0)
    m = rest.shape[0]
    cost = np.zeros(length + 1)
    cost[-1] = 1.0
    ones = np.ones((m, 1))
    a_ub = np.vstack([np.hstack([rest, -ones]), np.hstack([-rest, -ones])])
    b_ub = np.zeros(2 * m)
    a_eq = np.hstack([big[c], [0.0]])[None, :]
    bounds = [(None, None)] * length + [(0, None)]
    res = optimize.linprog(cost, a_ub, b_ub, a_eq, [1.0], bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    return _to_reference(res.x[:length].astype(complex))
