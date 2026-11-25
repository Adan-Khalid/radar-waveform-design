"""Pulse waveforms at a given time-bandwidth product, as complex baseband samples.

Conventions (used by every module):
- A pulse lasts T seconds and occupies a nominal bandwidth B, so its time-bandwidth
  product is BT. It is sampled at `samples_per_cell` samples per 1/B, which gives
  round(BT * samples_per_cell) samples.
- Every waveform is scaled to unit energy: sum(|s|^2) = 1.
- Phase codes have N chips of length T/N, so B = N/T and BT = N.
"""

from __future__ import annotations

import numpy as np

BARKER = {
    2: [1, -1],
    3: [1, 1, -1],
    4: [1, 1, -1, 1],
    5: [1, 1, 1, -1, 1],
    7: [1, 1, 1, -1, -1, 1, -1],
    11: [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1],
    13: [1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1],
}


def unit_energy(s: np.ndarray) -> np.ndarray:
    """Scale a waveform so that sum(|s|^2) = 1."""
    s = np.asarray(s, dtype=complex)
    return s / np.sqrt(np.sum(np.abs(s) ** 2))


def _sample_times(n: int) -> np.ndarray:
    """Sample instants as fractions of the pulse length, at the centre of each sample."""
    return (np.arange(n) + 0.5) / n


def lfm(time_bandwidth: float, samples_per_cell: int = 1) -> np.ndarray:
    """Linear FM pulse sweeping -B/2 to +B/2: s(t) = exp(j pi (B/T) (t - T/2)^2)."""
    n = int(round(time_bandwidth * samples_per_cell))
    t = _sample_times(n)  # in units of T
    return unit_energy(np.exp(1j * np.pi * time_bandwidth * (t - 0.5) ** 2))


def barker_phases(length: int = 13) -> np.ndarray:
    """Barker code as phases (0 or pi)."""
    if length not in BARKER:
        raise ValueError(f"no Barker code of length {length}; lengths: {sorted(BARKER)}")
    return np.where(np.array(BARKER[length]) > 0, 0.0, np.pi)


def frank_phases(length: int) -> np.ndarray:
    """Frank code of length M^2: phase(i, k) = 2 pi i k / M, chip index n = i M + k."""
    m = int(round(np.sqrt(length)))
    if m * m != length:
        raise ValueError("a Frank code needs a square length (4, 9, 16, ..., 64, ...)")
    i, k = np.meshgrid(np.arange(m), np.arange(m), indexing="ij")
    return (2 * np.pi * i * k / m).ravel()


def p4_phases(length: int) -> np.ndarray:
    """P4 code (Lewis and Kretschmer): phase_n = pi n^2 / N - pi n, n = 0..N-1."""
    n = np.arange(length)
    return np.pi * n**2 / length - np.pi * n


def phase_code(phases: np.ndarray, samples_per_chip: int = 1) -> np.ndarray:
    """Rectangular chips with the given phases, `samples_per_chip` samples each."""
    chips = np.exp(1j * np.asarray(phases, dtype=float))
    return unit_energy(np.repeat(chips, samples_per_chip))


def qpsk_symbols(n: int, seed: int) -> np.ndarray:
    """n random QPSK symbols exp(j(pi/4 + k pi/2)) from a seeded generator."""
    rng = np.random.default_rng(seed)
    return np.exp(1j * (np.pi / 4 + np.pi / 2 * rng.integers(0, 4, n)))


def ofdm(symbols: np.ndarray, samples_per_cell: int = 1) -> np.ndarray:
    """OFDM pulse: N subcarriers spaced 1/T, centred on zero, one symbol each, no cyclic prefix.

    s(t) = sum_k X_k exp(j 2 pi (k - (N-1)/2) t / T), 0 <= t < T, so B = N/T and BT = N.
    """
    symbols = np.asarray(symbols, dtype=complex)
    n_sub = symbols.size
    t = _sample_times(n_sub * samples_per_cell)
    freq = np.arange(n_sub) - (n_sub - 1) / 2  # in units of 1/T
    return unit_energy(np.exp(2j * np.pi * np.outer(t, freq)) @ symbols)


def papr_db(s: np.ndarray) -> float:
    """Peak-to-average power ratio of the samples, dB."""
    p = np.abs(np.asarray(s)) ** 2
    return float(10 * np.log10(p.max() / p.mean()))


def make_waveform(spec: dict, time_bandwidth: int, samples_per_cell: int = 1) -> np.ndarray:
    """Build a waveform from a settings entry such as {name: p4, length: 64}.

    The code length (or number of subcarriers) must equal the time-bandwidth product,
    so every waveform of one design point has the same B and T.
    """
    name = spec["name"]
    if name == "lfm":
        return lfm(time_bandwidth, samples_per_cell)
    if name == "ofdm":
        n = spec.get("subcarriers", time_bandwidth)
        _check_length(name, n, time_bandwidth)
        if spec.get("modulation", "qpsk") != "qpsk":
            raise ValueError("only QPSK subcarriers are implemented")
        return ofdm(qpsk_symbols(n, spec.get("seed", 2026)), samples_per_cell)
    phases = {"barker": barker_phases, "frank": frank_phases, "p4": p4_phases}
    if name not in phases:
        raise ValueError(f"unknown waveform {name!r}")
    n = spec.get("length", time_bandwidth)
    _check_length(name, n, time_bandwidth)
    return phase_code(phases[name](n), samples_per_cell)


def label(spec: dict, time_bandwidth: int) -> str:
    """Short display name, e.g. 'P4 64' or 'LFM (BT 64)'."""
    name = spec["name"]
    if name == "lfm":
        return f"LFM (BT {time_bandwidth})"
    n = spec.get("length", spec.get("subcarriers", time_bandwidth))
    pretty = {"barker": "Barker", "frank": "Frank", "p4": "P4", "ofdm": "OFDM"}[name]
    return f"{pretty} {n}"


def _check_length(name, n, time_bandwidth):
    if n != time_bandwidth:
        raise ValueError(
            f"{name}: length {n} differs from the time-bandwidth product {time_bandwidth}"
        )
