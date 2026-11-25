import numpy as np
import pytest

from radarwave import waveforms as wf
from radarwave.matched import autocorrelation, correlate
from radarwave.metrics import isl_db, mismatch_loss_db, psl_db, pulse_metrics, width


@pytest.mark.parametrize(
    "s",
    [
        wf.lfm(64, 8),
        wf.phase_code(wf.barker_phases(13), 8),
        wf.phase_code(wf.frank_phases(64), 1),
        wf.phase_code(wf.p4_phases(13), 4),
        wf.ofdm(wf.qpsk_symbols(64, 1), 2),
    ],
)
def test_unit_energy(s):
    assert np.isclose(np.sum(np.abs(s) ** 2), 1.0)


def test_lengths_and_constant_envelope():
    assert wf.lfm(13, 8).size == 104
    assert wf.phase_code(wf.p4_phases(64), 8).size == 512
    for s in (wf.lfm(64, 8), wf.phase_code(wf.p4_phases(64), 8)):
        assert np.allclose(np.abs(s), np.abs(s[0]))
    assert wf.papr_db(wf.lfm(64, 1)) == pytest.approx(0.0, abs=1e-12)
    assert wf.papr_db(wf.ofdm(wf.qpsk_symbols(64, 1), 1)) > 3.0


def test_frank_and_p4_phases_match_the_formulas():
    m = 8
    ph = wf.frank_phases(m * m)
    assert ph[3 * m + 5] == pytest.approx(2 * np.pi * 3 * 5 / m)
    n = 13
    p4 = wf.p4_phases(n)
    assert p4[7] == pytest.approx(np.pi * 49 / 13 - 7 * np.pi)
    with pytest.raises(ValueError):
        wf.frank_phases(13)


def test_lfm_instantaneous_frequency_spans_the_band():
    bt, os_ = 64, 8
    s = wf.lfm(bt, os_)
    f = np.diff(np.unwrap(np.angle(s))) / (2 * np.pi)  # cycles per sample
    f_cells = f * os_  # in units of B
    assert f_cells[0] == pytest.approx(-0.5, abs=0.01)
    assert f_cells[-1] == pytest.approx(0.5, abs=0.01)


def test_ofdm_energy_sits_on_its_subcarriers():
    n = 16
    sym = wf.qpsk_symbols(n, 3)
    s = wf.ofdm(sym, 1)
    # remove the half-subcarrier offset and half-sample timing, then a DFT gives the symbols
    t = np.arange(n)
    shift = np.exp(-2j * np.pi * (-(n - 1) / 2) * (t + 0.5) / n)
    spec = np.fft.fft(s * shift) * np.exp(-1j * np.pi * np.arange(n) / n)
    assert np.allclose(spec / spec[0], sym / sym[0], atol=1e-10)


def test_make_waveform_checks_lengths():
    s = wf.make_waveform({"name": "frank", "length": 64}, 64, 1)
    assert s.size == 64
    with pytest.raises(ValueError):
        wf.make_waveform({"name": "p4", "length": 13}, 64, 1)


def test_echo_peak_at_its_delay():
    s = wf.phase_code(wf.p4_phases(32), 1)
    x = np.zeros(200, complex)
    x[57 : 57 + 32] = s
    y, lag = correlate(x, s)
    assert lag[np.argmax(np.abs(y))] == 57
    yv, lagv = correlate(x, s, "valid")
    assert lagv[np.argmax(np.abs(yv))] == 57
    assert np.allclose(yv, y[lag >= 0][: yv.size])


def test_barker13_sidelobes_exact():
    for os_ in (1, 8):
        s = wf.phase_code(wf.barker_phases(13), os_)
        y, _ = autocorrelation(s)
        assert psl_db(y) == pytest.approx(20 * np.log10(1 / 13), abs=1e-9)
        assert isl_db(y) == pytest.approx(10 * np.log10(12 / 169), abs=1e-9)


def test_lfm_large_bt_approaches_sinc():
    y, _ = autocorrelation(wf.lfm(1000, 8))
    assert psl_db(y) == pytest.approx(-13.26, abs=0.3)
    assert isl_db(y) == pytest.approx(-9.68, abs=0.1)
    assert width(y, -3, 8) == pytest.approx(0.886, abs=0.01)


def test_rectangular_chip_mainlobe_width():
    # triangular mainlobe 1 - |tau|/tc: full width at -3 dB is 2 (1 - 10^(-3/20)) = 0.584 chip
    y, _ = autocorrelation(wf.phase_code(wf.barker_phases(13), 64))
    assert width(y, -3, 64) == pytest.approx(2 * (1 - 10 ** (-3 / 20)), abs=1e-9)


def test_mismatch_loss():
    s = wf.lfm(64, 4)
    assert mismatch_loss_db(s, s) == pytest.approx(0.0, abs=1e-12)
    assert mismatch_loss_db(s, s * np.hanning(s.size)) > 1.0
    m = pulse_metrics(s, samples_per_cell=4)
    assert m["mismatch_loss_db"] == pytest.approx(0.0, abs=1e-12)


def test_snr_gain_equals_time_bandwidth():
    rng = np.random.default_rng(5)
    for s in (wf.lfm(64, 1), wf.phase_code(wf.p4_phases(64), 1), wf.ofdm(wf.qpsk_symbols(64, 2))):
        n = s.size
        a = 0.3  # echo amplitude; per-sample SNR = a^2 / n with unit noise power
        w = (rng.standard_normal((20000, n)) + 1j * rng.standard_normal((20000, n))) / np.sqrt(2)
        y_noise = w @ np.conj(s)
        out_snr = a**2 / np.var(y_noise)
        gain_db = 10 * np.log10(out_snr / (a**2 / n))
        assert gain_db == pytest.approx(10 * np.log10(n), abs=0.1)
