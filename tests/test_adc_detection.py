import numpy as np
import pytest

from radarwave import theory
from radarwave import waveforms as wf
from radarwave.adc import (
    gaussian_optimal_full_scale,
    quantize,
    quantize_real,
    quantize_with_dither,
    staircase,
    step,
)
from radarwave.detection import (
    Receiver,
    TrialSet,
    crossing_db,
    detections,
    noise_statistics,
    run_receivers,
    thresholds_for,
    wilson,
)

# ---------------------------------------------------------------- quantizer


@pytest.mark.parametrize("bits", [1, 2, 3, 4])
def test_midrise_levels_symmetry_and_clipping(bits):
    a = 1.0
    x = np.linspace(-3, 3, 100001)
    q = quantize_real(x, bits, a)
    levels = np.unique(q)
    assert levels.size == 2**bits
    assert np.allclose(levels, -levels[::-1])
    assert levels.max() == pytest.approx(a - step(bits, a) / 2)
    thr, lev = staircase(bits, a)
    assert np.allclose(lev, levels)
    assert thr.size == 2**bits - 1


def test_one_bit_is_the_sign():
    x = np.array([-2.0, -1e-9, 1e-9, 5.0])
    assert np.array_equal(np.sign(quantize_real(x, 1, 3.0)), np.sign(x))


def test_midtread_has_a_zero_level():
    q = quantize_real(np.linspace(-2, 2, 10001), 2, 1.0, "midtread")
    assert np.allclose(np.unique(q), [-0.5, 0.0, 0.5])
    with pytest.raises(ValueError):
        quantize_real(np.zeros(3), 1, 1.0, "midtread")


def test_nonsubtractive_dither_makes_the_mean_linear():
    rng = np.random.default_rng(1)
    bits, a = 2, 1.0
    for x0 in (-0.6, -0.1, 0.03, 0.4):  # all inside |x| <= A - step/2 = 0.75
        u = rng.random(400_000) - 0.5
        z = quantize_with_dither(np.full(u.size, x0), bits, a, "midrise", "nonsubtractive", u)
        if abs(x0) <= a - step(bits, a) / 2:
            assert z.mean() == pytest.approx(x0, abs=3e-3)


def test_subtractive_dither_error_is_uniform_and_independent():
    rng = np.random.default_rng(2)
    bits, a = 3, 1.0
    x = rng.uniform(-0.8, 0.8, 400_000)
    u = rng.random(x.size) - 0.5
    e = quantize_with_dither(x, bits, a, "midrise", "subtractive", u) - x
    d = step(bits, a)
    assert np.var(e) == pytest.approx(d**2 / 12, rel=0.02)
    assert abs(np.corrcoef(e, x)[0, 1]) < 0.01


def test_complex_quantize_works_per_component():
    z = np.array([0.3 - 0.7j, -1.2 + 0.1j])
    q = quantize(z, 1, 1.0)
    assert np.allclose(q, 0.5 * (np.sign(z.real) + 1j * np.sign(z.imag)))
    assert np.array_equal(quantize(z, None, 1.0), z)


# ---------------------------------------------------------------- statistics helpers


def test_wilson_interval():
    lo, hi = wilson(90, 100)
    assert lo == pytest.approx(0.8256, abs=1e-3)
    assert hi == pytest.approx(0.9448, abs=1e-3)


def test_crossing_interpolates_and_flags():
    snr = np.array([0.0, 1.0, 2.0])
    assert crossing_db(snr, np.array([0.5, 0.8, 1.0]), 0.9) == pytest.approx(1.5)
    assert crossing_db(snr, np.array([0.1, 0.2, 0.3]), 0.9) == np.inf
    assert crossing_db(snr, np.array([0.95, 0.97, 1.0]), 0.9) == -np.inf


# ---------------------------------------------------------------- detection


def test_unquantized_threshold_and_pd_match_marcum():
    s = wf.phase_code(wf.p4_phases(64))
    pfa = 1e-3
    rx = Receiver()
    stats = noise_statistics(s, [rx], 400_000, seed=11)[rx]
    # noise-only |y|^2 is exponential with mean 1: threshold -ln(Pfa)
    thr = thresholds_for({rx: stats}, pfa)[rx]
    assert thr.value == pytest.approx(-np.log(pfa), rel=0.05)
    assert 0.0 <= thr.gamma <= 1.0
    snr = theory.required_snr_db(0.9, pfa)
    hits = detections(s, rx, thr, [snr], TrialSet(64, 6000, seed=12))
    lo, hi = wilson(hits.sum(), hits.size)
    assert lo - 0.01 <= 0.9 <= hi + 0.01


def test_one_bit_loss_near_two_db():
    s = wf.lfm(128, 1)
    rows, curves_rows, pfa_rows, _ = run_receivers(
        s,
        [Receiver(), Receiver(bits=1)],
        pfa=1e-3,
        pd_target=0.9,
        noise_stats=200_000,
        trials=4000,
        coarse_trials=1000,
        bootstrap=200,
        seed=3,
    )
    unq, one = rows
    assert unq["loss_db"] == 0.0
    assert unq["required_snr_db"] == pytest.approx(10.76, abs=0.4)
    assert one["loss_db"] == pytest.approx(2.0, abs=0.6)
    assert one["loss_lo"] < one["loss_db"] < one["loss_hi"]
    assert all(r["pfa_measured"] < 3e-3 for r in pfa_rows)
    assert len(curves_rows) > 10


def test_randomised_threshold_hits_pfa_for_one_bit_barker():
    # one bit and a binary code: the statistic takes few values, so ties need randomising
    s = wf.phase_code(wf.barker_phases(13))
    rx = Receiver(bits=1)
    pfa = 1e-3
    stats = noise_statistics(s, [rx], 300_000, seed=21)[rx]
    assert np.unique(np.round(stats, 3)).size < 100
    thr = thresholds_for({rx: stats}, pfa)[rx]
    assert 0.0 < thr.gamma < 1.0
    fresh = noise_statistics(s, [rx], 300_000, seed=22)[rx]
    u = np.random.default_rng(0).random(fresh.size)
    rate = thr.decide(fresh, u).mean()
    lo, hi = wilson(rate * fresh.size, fresh.size)
    assert lo - 2e-4 <= pfa <= hi + 2e-4


def test_receiver_labels_and_theory():
    assert Receiver().label == "unquantized"
    rx = Receiver(bits=3, dither="subtractive")
    assert "3 bit" in rx.label and "subtractive" in rx.label
    assert rx.a_over_sigma == pytest.approx(gaussian_optimal_full_scale(3))
    assert Receiver(bits=1).theory_loss_db() == pytest.approx(1.9612, abs=1e-3)
