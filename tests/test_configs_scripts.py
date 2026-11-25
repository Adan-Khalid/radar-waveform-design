import importlib.util
import numbers
from pathlib import Path

import pytest

from radarwave.config import REPO_ROOT, load_config, read_csv

SCRIPTS = REPO_ROOT / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_numbers_in_configs_are_numbers():
    # PyYAML reads 1.0e6 (no sign in the exponent) as text; every number must be a number.
    wf = load_config("waveforms")
    assert isinstance(wf["example_radar"]["carrier_hz"], numbers.Real)
    assert isinstance(wf["example_radar"]["bandwidth_hz"], numbers.Real)
    adc = load_config("adc")
    for key in ("pfa", "pd_target", "noise_statistics", "trials_per_snr", "fine_step_db"):
        assert isinstance(adc["common"][key], numbers.Real), key
    assert isinstance(load_config("sidelobes")["mismatched_filter"]["diagonal_loading"], float)


def test_design_points_share_bandwidth_and_duration():
    wf = load_config("waveforms")
    for dp in wf["design_points"].values():
        bt = dp["time_bandwidth"]
        for spec in dp["waveforms"]:
            assert spec.get("length", spec.get("subcarriers", bt)) == bt


@pytest.mark.parametrize("study", ["a1", "a2", "a3"])
def test_adc_study_quick_run(study, tmp_path: Path):
    _load("run_adc").main(["--study", study, "--quick", "--out", str(tmp_path)])
    rows = read_csv(tmp_path / f"{study}_losses.csv")
    assert rows and rows[0]["receiver"] == "unquantized"
    assert float(rows[0]["loss_db"]) == 0.0
    assert (tmp_path / f"{study}_pd_curves.csv").exists()
    pfa = read_csv(tmp_path / f"{study}_pfa_check.csv")
    assert all(float(r["pfa_measured"]) < 5e-3 for r in pfa)


@pytest.mark.parametrize("study", ["pd", "a1rep", "phase"])
def test_other_studies_quick_run(study, tmp_path: Path):
    _load("run_adc").main(["--study", study, "--quick", "--out", str(tmp_path)])
    curves = read_csv(tmp_path / f"{study}_pd_curves.csv")
    assert curves
    if study == "a1rep":
        assert read_csv(tmp_path / "a1rep_losses.csv")[0]["receiver"] == "unquantized"


def test_phase_study_shows_the_one_bit_barker_blind_spot():
    # One bit, Barker 13, Pfa 1e-4: an echo at phase 0 puts everything into I and is
    # rarely detected however strong it is; at 45 degrees it is found. P4 is found at both.
    # (At Pfa 1e-3 the threshold is low enough that the phase-0 echo gets through.)
    mod = _load("run_adc")
    cfg = load_config("adc")
    opts = mod.settings(cfg, quick=False)
    opts.update(noise_stats=300_000, trials=2000)
    _, rows, pfa = mod.study_phase(cfg, opts, quick=True)
    rows = [{k: str(v) for k, v in r.items()} for r in rows]

    def pd(wave, phase, snr):
        r = [
            x
            for x in rows
            if x["waveform"] == wave
            and float(x["echo_phase_deg"]) == phase
            and float(x["snr_db"]) == snr
        ]
        return float(r[0]["pd"])

    assert pd("Barker 13", 0.0, 25.0) < 0.3
    assert pd("Barker 13", 45.0, 25.0) > 0.95
    assert pd("P4 13", 0.0, 25.0) > 0.95
    assert pd("P4 13", 45.0, 25.0) > 0.95
    assert all(float(r["pfa_measured"]) < 3e-4 for r in pfa)
