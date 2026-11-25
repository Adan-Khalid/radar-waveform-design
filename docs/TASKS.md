# Tasks

Status: todo, done (written and run), checked (results reviewed against theory).

## Setup

| ID | Task | Done when | Status |
| --- | --- | --- | --- |
| T01 | Repository layout, pyproject, CI file, licence, citation file | folders in place, placeholder test passes | done |
| T02 | Research log, data card, plan, decisions | files written | done |
| T03 | Configs: `waveforms.yaml`, `sidelobes.yaml`, `adc.yaml`; `radarwave.config` | configs load in a test | done |
| T04 | Plot style `assets/style.mplstyle`, readable on light and dark backgrounds | one test figure looks right on both | done |

## Theory first

| ID | Task | Done when | Status |
| --- | --- | --- | --- |
| T10 | `radarwave.theory`: Marcum Q Pd, required SNR, window loss, quantizer efficiency (Bussgang) for bits, full scale, mid-rise or mid-tread, dither type | tests: 11.75 dB at Pd 0.9 and Pfa 1e-4; 1-bit 2/π; uniform 2, 3, 4 bits 0.55, 0.17, 0.05 dB at Max's step sizes | done |
| T11 | `LEARNING.md`: matched filter and SNR gain; PSL and ISL; ambiguity function; windows; mismatched filter; quantization and Bussgang; dither; detection with intervals | each section has maths, design choice, rejected alternatives, three interview questions | done |

## Waveforms and pulse compression

| ID | Task | Done when | Status |
| --- | --- | --- | --- |
| T20 | `radarwave.waveforms`: LFM, Barker, Frank, P4, OFDM, unit energy, samples per 1/B | tests: unit energy, constant envelope where expected, Frank and P4 phases against the formulas, OFDM spectrum on its subcarriers | done |
| T21 | `radarwave.matched` and `radarwave.metrics` | tests: Barker PSL −22.28 dB and ISL −11.49 dB, LFM at BT 1000 PSL near −13.26 dB, peak at the echo delay, mismatch loss 0 dB for the matched filter | done |
| T22 | SNR gain by Monte Carlo | measured gain within 0.1 dB of 10 log10(BT) for each waveform | done |
| T23 | `scripts/run_waveforms.py`: metrics table for both design points, OFDM spread over 200 draws | `reports/waveforms/metrics.csv`, `ofdm_spread.csv`, figures | done |

## Ambiguity and Doppler

| ID | Task | Done when | Status |
| --- | --- | --- | --- |
| T30 | `radarwave.ambiguity` (auto and cross) | tests: χ(0,0), zero-Doppler cut = matched filter, sinc zero-delay cut, discrete volume identity, LFM ridge slope and height | done |
| T31 | Doppler tolerance numbers and figures | `reports/ambiguity/doppler_tolerance.csv`, ambiguity and tolerance figures | done |

## Sidelobe control

| ID | Task | Done when | Status |
| --- | --- | --- | --- |
| T40 | Windows (own Taylor formula, checked against SciPy) on LFM | tests: window SNR loss equals (mean w)²/mean(w²); Hamming 1.34 dB, Hann 1.76 dB | done |
| T41 | Least-squares mismatched filter | tests: peak constraint met, ISL lower than matched, loss computed both ways agrees | done |
| T42 | `scripts/run_sidelobes.py`: tables and figures, filter length sweep, Doppler sensitivity | `reports/sidelobes/*.csv`, figures | done |
| T43 | Optional: minimum-PSL filter for Barker 13 by linear programming | compared in the table | done |

## Low-bit ADC study

| ID | Task | Done when | Status |
| --- | --- | --- | --- |
| T50 | `radarwave.adc`: mid-rise and mid-tread quantizers, full-scale rules, dither (none, non-subtractive, subtractive) | tests: level count, symmetry, clipping, one bit equals sign, dithered mean equals input inside the range | done |
| T51 | `radarwave.detection`: threshold from noise, Pd with Wilson intervals, coarse-to-fine required SNR, paired bootstrap | tests: Pfa on fresh noise, unquantized Pd inside the Marcum Q band | done |
| T52 | Run A1 (bits, dither, waveforms, BT 13 and 64) | `reports/adc/a1_*.csv`, Pd curves and loss figure | done |
| T53 | Run A2 (full scale, mid-rise and mid-tread, dither) | `reports/adc/a2_*.csv`, figure with theory lines | done |
| T54 | Run A3 (BT 13 to 1024) | `reports/adc/a3_*.csv`, figure | done |
| T55 | Optional: run A4 (weak target under a strong echo) | `reports/adc/a4_*.csv` | todo (not started; next step) |
| T56 | Replicate of A1 for LFM BT 64 with a new seed | `reports/adc/a1rep_*.csv` | done |
| T57 | Wide Pd curves against Marcum Q for the figure | `reports/adc/pd_pd_curves.csv` | done |

## MATLAB model

| ID | Task | Done when | Status |
| --- | --- | --- | --- |
| T60 | `matlab/+radarwave` functions mirroring the Python package | written | done (written, not run: needs MATLAB Online) |
| T61 | Unit tests and studies (waveforms, ambiguity, sidelobes, ADC with fewer trials), `run_all.m` | written | done (written, not run: needs MATLAB Online) |
| T62 | Python export of reference numbers for the cross-check | `reports/crosscheck/*.csv` | done |
| T63 | Run `run_all` in MATLAB Online | tests pass, cross-check PASS, output saved to `reports/matlab/` | todo (needs MATLAB Online) |

## Write-up

| ID | Task | Done when | Status |
| --- | --- | --- | --- |
| T70 | README with traced numbers, Mermaid diagram, figures with alt text and captions | every number checked against `reports/` | done (tables checked by script against the CSVs) |
| T71 | Hero figure and 1280 x 640 social preview | `assets/hero.png`, `assets/social-preview.png` | done (`scripts/make_readme_assets.py`) |
| T72 | Fresh-copy Quickstart test | commands run in a clean copy | done (2026-11-25: uv sync 95 s, 62 tests, ruff clean; the Quickstart scripts reproduced all 36 report CSVs byte for byte) |
