import numpy as np
import pytest
from scipy.signal import windows as sw

from radarwave import theory
from radarwave import waveforms as wf
from radarwave.adc import GAUSSIAN_OPTIMAL_STEP, gaussian_optimal_full_scale
from radarwave.ambiguity import ambiguity, doppler_tolerance
from radarwave.matched import correlate
from radarwave.metrics import isl_db, mismatch_loss_db, psl_db
from radarwave.sidelobes import (
    convolution_matrix,
    ls_mismatched_filter,
    minpsl_filter_real,
    taylor,
    window,
)

# ---------------------------------------------------------------- theory


def test_required_snr_marcum():
    assert theory.required_snr_db(0.9, 1e-4) == pytest.approx(11.749, abs=0.002)
    assert theory.pd_marcum(theory.required_snr_db(0.5, 1e-6), 1e-6) == pytest.approx(0.5)
    # Albersheim is an approximation of the same curve, within about 0.2 dB
    assert theory.albersheim_db(0.9, 1e-4) == pytest.approx(11.75, abs=0.2)


@pytest.mark.parametrize("bits,loss", [(1, 1.9612), (2, 0.5495), (3, 0.1657), (4, 0.0504)])
def test_quantizer_loss_at_max_steps(bits, loss):
    a = gaussian_optimal_full_scale(bits)
    assert theory.quantizer_loss_db(bits, a) == pytest.approx(loss, abs=2e-3)


def test_one_bit_is_two_over_pi_for_any_full_scale():
    for a in (0.1, 1.0, 30.0):
        assert theory.quantizer_efficiency(1, a) == pytest.approx(2 / np.pi, rel=1e-9)


def test_efficiency_against_brute_force_integration():
    x = np.linspace(-9, 9, 360001)
    pdf = np.exp(-(x**2) / 2) / np.sqrt(2 * np.pi)
    dx = x[1] - x[0]
    from radarwave.adc import quantize_real

    for bits, a, kind in [(2, 1.3, "midrise"), (3, 4.0, "midtread"), (4, 2.0, "midrise")]:
        q = quantize_real(x, bits, a, kind)
        brute = (np.sum(x * q * pdf) * dx) ** 2 / (np.sum(q * q * pdf) * dx)
        assert theory.quantizer_efficiency(bits, a, kind) == pytest.approx(brute, rel=1e-4)


def test_dither_efficiency_against_monte_carlo():
    rng = np.random.default_rng(0)
    x = rng.standard_normal(2_000_000)
    u = rng.random(x.size) - 0.5
    from radarwave.adc import quantize_with_dither

    for dither in ("nonsubtractive", "subtractive"):
        a = gaussian_optimal_full_scale(2)
        z = quantize_with_dither(x, 2, a, "midrise", dither, u)
        mc = np.mean(x * z) ** 2 / np.mean(z * z)
        assert theory.quantizer_efficiency(2, a, "midrise", dither) == pytest.approx(mc, rel=5e-3)


def test_max_steps_are_the_best_uniform_steps():
    for bits in (2, 3, 4):
        a0 = gaussian_optimal_full_scale(bits)
        best = theory.quantizer_efficiency(bits, a0)
        for f in (0.9, 1.1):
            assert theory.quantizer_efficiency(bits, a0 * f) < best
        assert GAUSSIAN_OPTIMAL_STEP[bits] > 0


def test_window_losses_harris():
    assert theory.window_loss_db(window("hamming", 4001)) == pytest.approx(1.34, abs=0.01)
    assert theory.window_loss_db(window("hann", 4001)) == pytest.approx(1.76, abs=0.01)


# ---------------------------------------------------------------- windows


def test_windows_match_scipy():
    for n in (13, 64, 512):
        assert np.allclose(window("hann", n), sw.hann(n, sym=True))
        assert np.allclose(window("hamming", n), sw.hamming(n, sym=True))
        assert np.allclose(window("blackman", n), sw.blackman(n, sym=True), atol=1e-15)
        assert np.allclose(taylor(n, 35, 5), sw.taylor(n, nbar=5, sll=35, norm=False))


# ---------------------------------------------------------------- ambiguity


def test_ambiguity_origin_and_zero_doppler_cut():
    s = wf.phase_code(wf.p4_phases(16), 2)
    chi, lags, nu = ambiguity(s, nu_t=[0.0, 0.5])
    assert abs(chi[lags == 0, 0][0]) == pytest.approx(1.0)
    y, ylags = correlate(s, s)
    assert np.allclose(chi[:, 0], y)
    assert np.array_equal(lags, ylags)


def test_zero_delay_cut_is_a_sampled_sinc():
    s = wf.lfm(32, 4)
    nu = np.linspace(-3, 3, 61)
    chi, lags, _ = ambiguity(s, nu_t=nu)
    expected = theory.dirichlet_zero_delay(nu, s.size)
    assert np.allclose(np.abs(chi[lags == 0][0]), expected, atol=1e-12)


def test_volume_identity_exact():
    s = wf.ofdm(wf.qpsk_symbols(13, 4), 2)
    ref = ls_mismatched_filter(wf.phase_code(wf.barker_phases(13)), 26)
    for r in (s, ref):
        k = 2 * s.size
        chi, _, _ = ambiguity(s, r, nu_t=np.arange(k) * s.size / k)
        vol = np.sum(np.abs(chi) ** 2) / k
        assert vol == pytest.approx(np.sum(np.abs(s) ** 2) * np.sum(np.abs(r) ** 2), rel=1e-10)


def test_lfm_ridge_height_and_slope():
    bt, os_ = 64, 8
    s = wf.lfm(bt, os_)
    for nut in (2.0, 5.0, 10.0):
        chi, lags, _ = ambiguity(s, nu_t=[nut])
        mag = np.abs(chi[:, 0])
        assert mag.max() == pytest.approx(theory.lfm_ridge_peak(nut, bt), abs=0.01)
        assert lags[np.argmax(mag)] / os_ == pytest.approx(-nut, abs=1 / os_ + 1e-9)


def test_doppler_tolerance_phase_code_vs_lfm():
    barker = wf.phase_code(wf.barker_phases(13), 4)
    summ, _ = doppler_tolerance(barker, nu_t_max=1.0, step=0.01, samples_per_cell=4)
    # constant-envelope pulse with no range-Doppler coupling: |sinc(nu T)| loses 3 dB near 0.44
    assert summ["nu_t_3db"] == pytest.approx(0.443, abs=0.02)
    assert summ["loss_at_0p5_db"] == pytest.approx(20 * np.log10(2 / np.pi), abs=0.05)
    lfm = wf.lfm(13, 4)
    summ_l, _ = doppler_tolerance(lfm, nu_t_max=1.0, step=0.01, samples_per_cell=4)
    assert summ_l["nu_t_3db"] == float("inf")
    assert summ_l["shift_at_1_cells"] == pytest.approx(-1.0, abs=0.25)


def test_p4_recovers_at_integer_doppler():
    p4 = wf.phase_code(wf.p4_phases(64), 4)
    summ, _ = doppler_tolerance(p4, nu_t_max=1.0, step=0.05, samples_per_cell=4)
    assert summ["loss_at_0p5_db"] < -3.0
    assert summ["loss_at_1_db"] > -0.5
    assert summ["shift_at_1_cells"] == pytest.approx(-1.0)


# ---------------------------------------------------------------- mismatched filters


def test_convolution_matrix():
    s = np.array([1.0, 2.0, 3.0])
    g = np.array([1.0, -1.0])
    assert np.allclose(convolution_matrix(s, 2) @ g, np.convolve(s, g))


def test_ls_filter_constraint_and_gain():
    code = wf.phase_code(wf.barker_phases(13))
    for k in (13, 39):
        ref = ls_mismatched_filter(code, k)
        y, lags = correlate(code, ref)
        c = (13 + k - 2) // 2
        assert y[lags == c - (k - 1)][0] == pytest.approx(1.0)
        assert isl_db(y) < 10 * np.log10(12 / 169)
        assert 0 < mismatch_loss_db(code, ref) < 1.0


def test_ls_filter_beats_matched_isl_for_p4():
    code = wf.phase_code(wf.p4_phases(64))
    y_mf, _ = correlate(code, code)
    y_ls, _ = correlate(code, ls_mismatched_filter(code, 192))
    assert isl_db(y_ls) < isl_db(y_mf) - 5


def test_minpsl_filter_has_lower_peak_sidelobe_than_ls():
    code = wf.phase_code(wf.barker_phases(13))
    y_ls, _ = correlate(code, ls_mismatched_filter(code, 39))
    y_lp, _ = correlate(code, minpsl_filter_real(code, 39))
    assert psl_db(y_lp) <= psl_db(y_ls) + 1e-6
