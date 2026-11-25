# v1.0.0

First release: five radar pulses compared at equal bandwidth and duration, and a measured answer to how much detection a receiver gives up when its ADC keeps only 1 to 4 bits per sample, checked against theory.

## What it does

- LFM, Barker 13, Frank, P4 and OFDM pulses at the same bandwidth and duration (13 and 64 range cells), all at unit energy.
- Matched filtering with peak and integrated sidelobe levels, resolution and SNR gain; an own narrowband ambiguity function and Doppler tolerance for every pulse.
- Sidelobe control and its SNR cost: amplitude windows on LFM (Taylor written out from its formula), least-squares mismatched filters on phase codes and a minimum-peak-sidelobe filter for Barker 13 by linear programming.
- A 1 to 4-bit I/Q ADC model (mid-rise and mid-tread, with and without dither) and a Monte Carlo detector at Pfa 10⁻⁴ with a randomised threshold, common random numbers, Wilson intervals and a paired bootstrap, compared with the Bussgang efficiency of each quantizer.
- 62 tests against closed-form theory and small Monte Carlo runs, run on every push by GitHub Actions.

## Key results

Pd 0.9 at Pfa 10⁻⁴, simulation (`reports/summary/headline.csv` and the files it cites):

| Measure | Reference | This release |
| --- | --- | --- |
| SNR needed with no quantization (LFM, time-bandwidth 64) | 11.75 dB (Marcum Q) | 11.75 dB [11.72, 11.79] |
| Loss of a 3-bit ADC (same pulse) | 0.17 dB (theory) | 0.17 dB [0.15, 0.19] |
| Loss of a 1-bit ADC (same pulse) | 1.96 dB (theory) | 2.34 dB [2.30, 2.38] |
| Loss of a 1-bit ADC, P4 code of 1024 and 13 chips | 1.96 dB (theory) | 2.00 dB and 5.43 dB |
| Weak target when a strong echo sets the gain (3 bits, full scale 30 dB above the noise) | plain mid-rise ADC: 2.34 dB | dither removed: 7.90 dB; mid-tread: 11.52 dB |
| Barker 13 peak sidelobe | −22.3 dB (matched filter) | −38.5 dB with a 39-tap least-squares filter, for 0.20 dB |

A full rerun with new random numbers moves the losses by up to 0.08 dB (`reports/adc/a1_replicate_check.csv`).

## How to run

```bash
uv sync --group dev
uv run pytest
uv run python scripts/run_adc.py --study a1
```

The README lists every script, what it produces and how long it takes.

## Known limits

- Simulation only: ideal ADC apart from quantization and clipping, white Gaussian noise, non-moving point targets with constant amplitude, one complex sample per range cell, a detector that knows the target's cell.
- The strong-echo study assumes the strong echo does not overlap the weak target's samples.
- The MATLAB version (`matlab/`) is written but has not been run yet; every number comes from the Python code.
