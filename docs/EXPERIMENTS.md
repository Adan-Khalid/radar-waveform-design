# Experiments

One entry per run that produces a number used anywhere: script, settings, seed, output file and result. Theory first, then the waveform studies, then the ADC study (baseline first, one change at a time). Every run is deterministic: the scripts use fixed seeds (2026 unless stated), so running them again gives the same files.

## Theory (`scripts/run_theory.py`, no randomness)

| ID | What | Output | Result |
| --- | --- | --- | --- |
| T1 | Required SNR, one pulse, non-fluctuating target, Marcum Q and Albersheim | `reports/theory/required_snr.csv` | Pd 0.9 at Pfa 10⁻⁴: 11.749 dB (Albersheim 11.67 dB) |
| T2 | Window SNR loss, (mean w)² / mean(w²) | `reports/theory/window_loss.csv` | Hann 1.761, Hamming 1.344, Taylor (35 dB, n̄ 5) 0.926, Blackman 2.372 dB |
| T3 | Quantizer loss at Max's step, every bit depth, kind and dither | `reports/theory/quantizer_loss_optimal.csv` | mid-rise, no dither: 1.961, 0.550, 0.166, 0.050 dB for 1 to 4 bits; non-subtractive dither raises them to 2.84, 0.90, 0.29, 0.09 dB |
| T4 | Quantizer loss against full scale (curves for the A2 figure) | `reports/theory/quantizer_loss_vs_full_scale.csv` | 261 full-scale values from −15 to +50 dB, every bit depth, kind and dither |

## Waveforms (`scripts/run_waveforms.py`, `configs/waveforms.yaml`, 8 samples per cell)

| ID | What | Output | Result |
| --- | --- | --- | --- |
| W1 | PSL, ISL, −3 dB width, PAPR of every waveform | `reports/waveforms/metrics.csv` | BT 13: Barker −22.3 dB PSL, P4 −17.3, LFM −14.4, OFDM −11.3. BT 64: Frank −27.8, P4 −24.4, OFDM −13.7, LFM −13.5. Phase codes have the narrowest mainlobes (0.58 to 0.63 cells against 0.87 to 0.88 for LFM and OFDM) |
| W2 | Matched-filter SNR gain, Monte Carlo (20,000 trials, one sample per cell) | `reports/waveforms/snr_gain.csv` | within 0.08 dB of 10 log10(BT) for all eight waveforms (11.14 and 18.06 dB) |
| W3 | OFDM spread over 200 random QPSK draws | `reports/waveforms/ofdm_spread.csv` | median PSL −12.6 dB (13 subcarriers) and −13.6 dB (64); median ISL −4.3 and −3.5 dB; median PAPR of the samples taken at B: 5.0 and 6.8 dB (95th percentile 7.1 and 8.6 dB) |

## Doppler (`scripts/run_ambiguity.py`, νT from 0 to 4 in steps of 0.01, then to BT for LFM)

| ID | What | Output | Result |
| --- | --- | --- | --- |
| D1 | Doppler tolerance of every waveform | `reports/ambiguity/doppler_tolerance.csv`, `doppler_cuts.csv` | for every code the peak first drops 3 dB near νT = 0.44 (it follows the zero-delay cut, |sinc(νT)| for a constant envelope, until about νT = 0.5). P4 and Frank then recover at integer νT, one cell further on (P4 64 −0.14 dB at νT = 1, −0.62 dB at νT = 4; Frank 64 −0.36 and −1.64 dB). Barker 13 and OFDM do not recover (−7.0 and −12.2 dB at νT = 1 for Barker 13 and OFDM 64). LFM keeps its peak far longer (3 dB at νT 3.8 for BT 13 and 18.7 for BT 64, i.e. 0.29 B) but moves one range cell per unit of νT |

With the example radar (10 GHz, B = 1 MHz) one unit of νT is 234 m/s for the 64 µs pulse and 1,153 m/s for the 13 µs pulse (`speed_per_unit_nut_mps` in the same file).

## Sidelobe control (`scripts/run_sidelobes.py`, `configs/sidelobes.yaml`)

| ID | What | Output | Result |
| --- | --- | --- | --- |
| S1 | Receive windows on LFM, BT 13 and 64 (8 samples per cell) | `reports/sidelobes/lfm_windows.csv` | BT 64: Taylor −32.3 dB PSL for 0.93 dB, Hamming −38.6 dB for 1.35 dB, Hann −31.4 dB for 1.77 dB, Blackman −37.9 dB for 2.38 dB, against −13.5 dB with no window. BT 13: every window stops near −24 dB. The measured loss equals the formula in every case |
| S2 | Least-squares mismatched filters, lengths N to 4N, and a minimum-PSL filter for Barker 13 (one sample per chip) | `reports/sidelobes/mismatched_filters.csv` | Barker 13: PSL −22.3 dB (matched) to −38.5 dB with 3N taps for 0.20 dB and −48.4 dB with 4N taps for 0.21 dB; minimum-PSL filter with 3N taps −43.1 dB for 0.19 dB. Frank 64: −39.5 dB for 0.30 dB (3N). P4 64: −35.5 dB for 0.41 dB (3N). P4 13 gains least (−28.6 dB for 0.88 dB) |
| S3 | Sidelobes of the mismatched filters under Doppler | `reports/sidelobes/mismatched_doppler.csv` | the least-squares filters (3N taps) lose their advantage with Doppler: Barker 13 goes from −38.5 dB at νT = 0 to −32.7 dB at 0.05 and is no better than the matched filter from νT = 0.18; P4 13, Frank 64 and P4 64 from νT = 0.28 to 0.29 |

## Low-bit ADC (`scripts/run_adc.py`, `configs/adc.yaml`)

Common settings: one complex sample per range cell, Pfa 10⁻⁴, Pd target 0.9, threshold from 2 × 10⁶ noise-only statistics per receiver (randomised on ties, D9), false-alarm check on 10⁶ fresh statistics, coarse grid 0.5 dB with 4,000 trials, fine grid ±1 dB in 0.1 dB steps with 20,000 trials, paired Poisson bootstrap with 1,000 replicates, seed 2026. Loss = required SNR minus that of the unquantized receiver on the same trials. Settings of each run: `reports/adc/<run>_settings.csv`.

### E-A1a: first bit-depth run, plain threshold (superseded)

`reports/adc/superseded/a1_plain_threshold_*.csv`. Barker 13 at one bit needed 13.78 dB more SNR. The cause was first suspected to be the threshold: with one bit and a 13-chip binary code the statistic takes few values, so a plain threshold cannot give Pfa = 10⁻⁴ exactly. The randomised threshold (D9) was added and every run repeated (E-A1). The Barker 13 result did not change: it is a real effect of the echo phase (see E-A1, D12 and `LEARNING.md` section 8). All other numbers were identical, since their statistics are continuous.

### E-A1: bits and dither, every waveform (`--study a1`, 306 s)

Output: `reports/adc/a1_losses.csv`, `a1_pd_curves.csv`, `a1_pfa_check.csv`. Loss in dB against unquantized, 95 % interval in brackets for LFM; theory = Bussgang efficiency at Max's step.

| Receiver | LFM BT 64 | Frank 64 | P4 64 | OFDM 64 | LFM BT 13 | Barker 13 | P4 13 | OFDM 13 | Theory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| unquantized, required SNR | 11.75 | 11.75 | 11.77 | 11.77 | 11.74 | 11.78 | 11.77 | 11.75 | 11.75 |
| 4 bit | 0.07 [0.06, 0.09] | 0.07 | 0.05 | 0.05 | 0.16 | 0.13 | 0.12 | 0.20 | 0.05 |
| 3 bit | 0.17 [0.15, 0.19] | 0.20 | 0.19 | 0.18 | 0.40 | 0.38 | 0.36 | 0.59 | 0.17 |
| 2 bit | 0.60 [0.57, 0.63] | 0.66 | 0.67 | 0.72 | 1.20 | 1.17 | 1.14 | 1.67 | 0.55 |
| 1 bit | 2.34 [2.30, 2.38] | 2.38 | 2.42 | 2.78 | 5.52 | 13.78 | 5.43 | 9.97 | 1.96 |
| 3 bit, dither not removed | 0.33 [0.31, 0.35] | 0.33 | 0.30 | 0.31 | 0.55 | 0.52 | 0.46 | 0.67 | 0.29 |
| 3 bit, dither removed | 0.21 [0.19, 0.23] | 0.20 | 0.18 | 0.23 | 0.42 | 0.39 | 0.36 | 0.56 | 0.17 |
| 1 bit, dither not removed | 3.24 [3.20, 3.28] | 3.26 | 3.18 | 3.76 | 6.25 | 14.69 | 6.25 | 10.86 | 2.84 |
| 1 bit, dither removed | 3.11 [3.07, 3.15] | 3.07 | 3.05 | 3.78 | not reached | 18.09 | not reached | not reached | 2.55 |

Findings:
- At BT 64, without dither, the constant-envelope pulses (LFM, Frank, P4) lose within 0.12 dB of the theory at 2, 3 and 4 bits; one bit costs them 2.34 to 2.42 dB against 1.96 dB in theory, because each sample's SNR is not yet small (about −6.3 dB at the unquantized operating point, 11.75 − 18.06 dB, and −4.0 dB at the one-bit operating point).
- When noise sets the full scale, dither does not help, as the theory predicts: dither left in the samples costs 0.73 to 0.97 dB more at one bit, 0.33 to 0.51 dB at 2 bits, 0.09 to 0.16 dB at 3 bits and at most 0.06 dB at 4 bits; removing it after the ADC brings the loss back to within 0.06 dB of no dither at 3 and 4 bits, but not at 1 and 2 bits (0.06 to 0.35 dB more at 2 bits).
- The waveform hardly matters at BT 64 for constant-envelope pulses (spread under 0.1 dB from 2 to 4 bits). OFDM loses a little more at 1 and 2 bits (2.78 and 0.72 dB), as its peaks raise the SNR of some samples.
- At BT 13 everything costs more (3 bits about 0.4 dB, one bit about 5.5 dB). Barker 13 at one bit needs 13.8 dB more: a real code puts the whole echo into I or Q depending on its phase, and at phase 0 the one-bit statistic cannot reach the threshold. One bit with the dither removed never reaches Pd 0.9 at BT 13: at 60 dB Pd is 0.85 (LFM), 0.86 (P4) and 0.66 (OFDM) (`a1_pd_curves.csv`, D12).
- False-alarm check: measured Pfa between 0.66 × 10⁻⁴ and 1.24 × 10⁻⁴ for all 104 receivers; 84 of the 104 Wilson intervals contain 10⁻⁴. The rest reflect the error of the threshold itself (estimated from 2 × 10⁶ statistics), which that interval does not include. Because all receivers use the same noise, this error largely cancels in the losses; E-A1r checks it with a full replicate.

### E-PD: Pd curves for the figure (`--study pd`, 19 s)

`reports/adc/pd_pd_curves.csv`: LFM BT 64, no dither, 4 to 18 dB in 0.5 dB steps, 20,000 trials per point. The unquantized curve lies inside its 95 % Wilson interval around the Marcum Q prediction at 29 of 29 points.

### E-A3: code length (`--study a3`, 337 s)

`reports/adc/a3_losses.csv`: P4 code, no dither, Max's step, BT 13 to 1024.

| Bits | BT 13 | 32 | 64 | 128 | 256 | 1024 | Theory |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 5.43 | 2.80 | 2.42 | 2.21 | 2.12 | 2.00 [1.97, 2.04] | 1.96 |
| 2 | 1.14 | 0.71 | 0.67 | 0.60 | 0.61 | 0.54 [0.52, 0.57] | 0.55 |
| 3 | 0.36 | 0.21 | 0.19 | 0.17 | 0.16 | 0.16 [0.15, 0.18] | 0.17 |
| 4 | 0.12 | 0.05 | 0.05 | 0.02 | 0.07 | 0.07 [0.06, 0.09] | 0.05 |

The measured loss falls toward the low-SNR theory as the code gets longer and reaches it at BT 1024 (one bit: 2.00 dB, interval 1.97 to 2.04, against 1.96). The unquantized required SNR is 11.71 to 11.78 dB at every length (theory 11.75).

### E-A1r: full replicate with new random numbers (`--study a1rep`, seed 2027, 67 s)

`reports/adc/a1rep_losses.csv`, comparison in `reports/adc/a1_replicate_check.csv` (from `scripts/summarise_adc.py`). LFM BT 64, all 13 receivers, new noise, phases, dither and thresholds. The losses moved by at most 0.08 dB (one bit with dither left in: 3.24 to 3.32 dB; 2 bits: 0.60 to 0.66 dB); 11 of 13 differences lie inside the two runs' combined intervals. The bootstrap intervals cover the trial-to-trial noise but not the error of the threshold itself, so the honest run-to-run uncertainty of a loss is about ±0.08 dB, slightly wider than the printed intervals (D13).

### E-A2: ADC full scale (`--study a2`)

`reports/adc/a2_losses.csv`, `a2_pd_curves.csv`, `a2_pfa_check.csv`; figure `reports/figures/adc_full_scale.png`. LFM BT 64; full scale A/σ per component from −12 to +46 dB in 2 dB steps (30 values); 1 to 4 bits; mid-rise and mid-tread; no dither, dither left in, dither removed; 21 receivers plus the unquantized reference at every full scale, 631 receivers in all. Run time about 50 minutes of computation; the settings file shows 80,331 s because the laptop slept during the run (two groups of receivers were interrupted overnight and finished on waking; sleep does not change the random numbers, so the results are unaffected).

Loss (dB) at selected full scales (the optimum for noise is about 4 to 9 dB):

| Receiver | −12 | 0 | 6 | 12 | 18 | 24 | 30 | 46 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 bit, no dither | 2.34 | 2.34 | 2.34 | 2.34 | 2.34 | 2.34 | 2.34 | 2.34 |
| 3 bit, mid-rise, no dither | 1.81 | 0.78 | 0.25 | 0.31 | 1.24 | 2.34 | 2.34 | 2.34 |
| 3 bit, mid-rise, dither removed | 1.86 | 0.78 | 0.26 | 0.32 | 1.25 | 3.59 | 7.90 | 23.09 |
| 3 bit, mid-rise, dither left in | 1.84 | 0.82 | 0.35 | 0.67 | 2.25 | 5.86 | 11.54 | 27.46 |
| 3 bit, mid-tread, no dither | 1.91 | 0.95 | 0.38 | 0.36 | 1.20 | 4.41 | 11.52 | 30.86 |
| 4 bit, mid-rise, no dither | 1.78 | 0.68 | 0.15 | 0.07 | 0.31 | 1.25 | 2.34 | 2.34 |

Findings:
- **A mid-rise quantizer without dither never loses more than its one-bit value**, however large the full scale: once the step is much larger than the noise, every sample falls in the two middle levels and the receiver becomes a one-bit (sign) receiver. All bit depths end at 2.34 dB, the one-bit loss of this pulse (A1).
- **Dither and mid-tread quantizers lose without limit** as the full scale grows, approaching 1 dB more loss for every dB of full scale once the step is well above the noise (dither adds power in proportion to the step squared; a mid-tread quantizer sends weak samples to zero). At 30 dB above the noise, 3 bits lose 7.9 dB with dither removed, 11.5 dB with dither left in and 11.5 dB with a mid-tread quantizer, against 2.3 dB for plain mid-rise. Dither does help a mid-tread quantizer (19.5 against 30.9 dB at 46 dB for 3 bits) but never beats plain mid-rise in this setting.
- Each extra bit moves the point where the loss starts to rise by 6 dB, as expected from halving the step.
- Theory (lines in the figure) and Monte Carlo agree wherever each sample's SNR stays low. Two visible gaps: one bit sits about 0.4 dB above its theory everywhere (finite per-sample SNR, as in A1), and mid-tread without dither at large full scale loses far less than the low-SNR theory, which goes to infinity; at those losses the echo needed for Pd 0.9 is no longer weak per sample and crosses the zero level by itself.
- False alarms: measured Pfa from 0.74 × 10⁻⁴ to 1.27 × 10⁻⁴ for all mid-rise and unquantized receivers (359 of 390 intervals contain 10⁻⁴). For mid-tread without dither at large full scale the noise output is almost always exactly zero, so the statistic is nearly discrete and the statistics of one noise record are strongly correlated; the achieved Pfa reaches 1.56 × 10⁻⁴. That makes those (5 to 31 dB) losses optimistic by at most about 0.17 dB (the change in required SNR between Pfa 10⁻⁴ and 1.56 × 10⁻⁴, `theory.required_snr_db`).
- The unquantized reference needs 11.75 dB at every full scale (same trials), as in A1.

What this does and does not cover: here the strong echo only sets the gain; it does not overlap the weak target's samples. When it does overlap, the quantization error becomes a deterministic function of the strong echo (spurs), and that is where dither is expected to help (open run A4).

### E-A2s: the full scale read as a strong echo (`scripts/summarise_adc.py`)

`reports/adc/a2_strong_echo.csv`. If the gain is set so that a strong echo's largest I or Q sample just reaches full scale, A/σ = √(2 · SNR_strong · PAPR), with SNR_strong its mean SNR per sample. The OFDM 64 pulse has a PAPR of 7.3 dB on its samples (`reports/waveforms/metrics.csv`), so it needs 7.3 dB more full scale for the same echo power. Loss interpolated from E-A2 (mid-rise, no dither):

| Strong echo SNR per sample | 10 dB | 20 dB | 30 dB |
| --- | --- | --- | --- |
| 3 bit, constant-envelope strong echo | 0.42 | 2.26 | 2.34 |
| 3 bit, OFDM strong echo | 1.80 | 2.34 | 2.34 |
| 4 bit, constant-envelope strong echo | 0.12 | 1.05 | 2.34 |
| 4 bit, OFDM strong echo | 0.57 | 2.34 | 2.34 |

An OFDM strong echo uses up its 7.3 dB of peak-to-average ratio in dynamic range: a 3-bit receiver loses 1.80 dB next to a 10 dB OFDM echo against 0.42 dB next to a constant-envelope echo of the same power.

### E-PH: one bit and the echo phase (`--study phase`, 3 s)

`reports/adc/phase_pd_curves.csv`: Barker 13 and P4 13 at one bit, the echo phase fixed at 0, 10, 22.5 and 45 degrees instead of random, SNR 12 to 40 dB, 20,000 trials per point, threshold from 2 × 10⁶ noise statistics (Pfa 10⁻⁴), same receiver and decision code as every other run (`detections` with a fixed phase). Barker 13 at phase 0 levels off at Pd 0.09 from 20 dB up to 40 dB, while at 45 degrees it is found with Pd 0.95 at 16 dB. P4 13 is found at every phase (Pd 0.69 to 0.88 at 16 dB, 0.97 to 1.00 at 20 dB). This confirms the explanation of the large one-bit loss of Barker 13 in E-A1. The ceiling depends on the false-alarm setting: at Pfa 10⁻³ the lower threshold lets the phase-0 echo through far more often.

Correction: a first version of this run (2026-11-23) quantized its trials with a different ADC scale from its threshold, which reported Pd 0.003 and 0.72 for Barker 13 and a blind spot of P4 at 45 degrees. Both were wrong; an independent review found the mismatch, the run now uses the shared detector, and a test checks it (`test_phase_study_shows_the_one_bit_barker_blind_spot`).
