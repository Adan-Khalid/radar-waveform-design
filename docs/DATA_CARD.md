# Data card

## Source

No external dataset. Every signal in this repository is simulated by the code in `src/radarwave/` from the settings in `configs/`. Nothing is downloaded.

| Item | Value |
| --- | --- |
| Signals | radar pulses (LFM, Barker 13, Frank, P4, OFDM), point echoes, white circular Gaussian noise |
| Generator | `numpy.random.default_rng` (PCG64) in Python; `rng(seed, 'twister')` in MATLAB |
| Seeds | stated in each config file and in `docs/EXPERIMENTS.md` for every run |
| Licence | the code and the generated numbers are MIT, like the rest of the repository |
| Size | the reported tables and figures are small (under 10 MB); Monte Carlo trials are generated on the fly and not stored |

## What the simulation assumes

- Complex baseband, narrowband: a Doppler shift is a phase ramp over the pulse; no pulse stretching.
- One point echo per target, constant amplitude over the pulse (non-fluctuating), random phase from trial to trial.
- Thermal noise white and circular Gaussian at the ADC sample rate. In the ADC study the sample rate equals the bandwidth, one complex sample per range cell.
- Ideal ADC apart from quantization and clipping: no jitter, no offset, no integral or differential nonlinearity, I and Q matched.
- Echoes from a strong object, where used, are point echoes too; no clutter distributions.

These are the textbook assumptions under which the theory in `docs/LEARNING.md` holds, so measured numbers can be checked against it. They are also the limits of the results: real receivers add phase noise, I/Q imbalance, ADC nonlinearity and clutter.

## MATLAB and Python random numbers

The two languages use different generators, so Monte Carlo numbers from MATLAB and Python are compared within their confidence intervals, not sample by sample. Deterministic quantities (waveforms, filters, sidelobe levels, ambiguity functions, theory values) must agree to 1e-9.
