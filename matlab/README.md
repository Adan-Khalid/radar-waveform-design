# MATLAB version

Status: written, not yet run (no MATLAB licence was available while it was written). Until it is run, every result in this repository comes from the Python code.

The `+radarwave` package implements the same equations as the Python package `src/radarwave`, with the same names in camelCase (`radarwave.lfm`, `radarwave.lsMismatchedFilter`, `radarwave.quantizerEfficiency`). It needs only base MATLAB, the Signal Processing Toolbox (`marcumq`, `taylorwin` in the tests), the Statistics and Machine Learning Toolbox (`poissrnd`, `prctile`, `ncx2cdf`, `corr`) and the Optimization Toolbox (`linprog`, for the minimum-PSL filter). All four are in the free MATLAB Online licence.

| Folder | What it holds |
| --- | --- |
| `+radarwave/` | waveforms, matched filter, metrics, ambiguity function, windows, mismatched filters, quantizer, theory, Monte Carlo detection |
| `tests/TestRadarwave.m` | unit tests against closed-form values (Barker sidelobes, sinc limits, volume identity, Taylor against `taylorwin`, quantizer losses, Marcum Q) |
| `studies/` | `study_waveforms`, `study_ambiguity`, `study_sidelobes`, `study_adc`, `study_crosscheck` |
| `run_all.m` | runs the tests, then every study |

## Running it in MATLAB Online

1. Upload the `matlab/` and `reports/` folders to MATLAB Drive, keeping the folder names (one parent folder, for example `radar-waveform-design/`). `reports/crosscheck/` holds the OFDM symbols, waveform samples and filter taps exported by `scripts/export_crosscheck.py`, because MATLAB and NumPy draw different random numbers.
2. In the command window: `cd radar-waveform-design/matlab`, then `run_all`.
3. The results go to `reports/matlab/`: CSV tables, PNG figures and `crosscheck.csv` with one PASS or FAIL line per comparison with Python.

`study_adc` runs the Monte Carlo detection study with fewer trials than Python (a few minutes). If a session stops it, lower `opts.noiseStats` and `opts.trials` in `studies/study_adc.m`, or run `run_all(true)` to skip it.
