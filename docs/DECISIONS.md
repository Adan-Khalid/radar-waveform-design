# Decisions

One entry per design decision: what was chosen, the alternatives, and why.

## D1 (2026-11-22): equal bandwidth and duration in every comparison

Waveforms are compared at the same bandwidth B and duration T, so at the same time-bandwidth product BT and the same nominal range resolution. Two design points: BT = 13, where Barker 13 exists, and BT = 64, where Frank (8 x 8) exists and the per-sample SNR is low enough for the low-bit theory to apply. P4 and OFDM exist at any length and appear in both.

"Bandwidth" here is the nominal one: the sweep width for LFM, 1/(chip length) for phase codes, N times the subcarrier spacing for OFDM. The spectra still differ in shape (phase codes have sinc-shaped spectra wider than B), which `LEARNING.md` explains; the alternative of matching a 3 dB or 99 % power bandwidth was rejected because it would give every code a different chip length and make the familiar textbook values (Barker −22.3 dB, LFM −13.3 dB) harder to check.

## D2 (2026-11-22): Frank and P4 both, not Frank or P4

P4 exists for any length, so it appears at BT = 13 and 64; Frank needs a square length and appears at 64. Having both at 64 shows the known difference in Doppler tolerance between them at no extra cost.

## D3 (2026-11-22): OFDM with random QPSK and no cyclic prefix

The OFDM pulse carries random QPSK symbols from a fixed seed, as a communication payload would, and has no cyclic prefix, so it is processed by the same matched filter as the other pulses. Its sidelobes depend on the symbols, so its metrics are also given as a spread over 200 symbol draws. Rejected for now: low-PAPR symbol choices (for example Newman phases) and the reciprocal (spectral division) receive filter used in OFDM radar; both are left for future work on OFDM and OTFS sensing.

## D4 (2026-11-22): one complex sample per range cell in the ADC study

The ADC study samples at f_s = B, so noise samples are independent and each sample is one range cell. This is where the quantizer theory is exact and the effect of the bit depth is clearest. Oversampled low-bit receivers (more samples per cell, which recover part of the loss) are listed as a next step. Waveform metrics and ambiguity functions use 8 samples per 1/B so sidelobe peaks between cells are not missed.

## D5 (2026-11-22): clairvoyant single-cell detector for the detection loss

The detection loss is measured with the threshold set from noise-only statistics of the same receiver, at the known target cell. This isolates the loss caused by the ADC from CFAR loss and straddle loss, which are separate topics. Pfa = 10⁻⁴, threshold from at least 10⁶ noise-only statistics, checked on fresh noise.

## D6 (2026-11-22): theory from the Bussgang efficiency, measurement from detection only

For every receiver the low-SNR theory gives an efficiency η (Bussgang), so the predicted loss is −10 log10 η, computed exactly for every bit depth, full scale and dither setting (`theory.quantizer_efficiency`). The measured result is always a detection run (Pd at fixed Pfa, required SNR for Pd 0.9), never a proxy. A cheaper effective-SNR (deflection) estimate was planned for the full-scale sweep, but simulating only the samples that feed the target cell made full detection runs fast enough (about 3 s per receiver), so the proxy was dropped. Because the matched-filter output sums 13 to 1024 samples it is close to Gaussian, so detection loss and SNR loss should agree when each sample's SNR is low; run A3 shows where that stops being true.

## D7 (2026-11-22): own implementations, written twice

Waveforms, matched filter, ambiguity function, windows (including the Taylor formula), mismatched filter, quantizer and detector are written from the equations in Python and in MATLAB, using base MATLAB with the Signal Processing and Statistics toolboxes (and Optimization for the optional linear-programming filter), so everything is meant to run in free MATLAB Online, though it has not been run yet (checked on the MathWorks page on 2026-11-22: it has those toolboxes but not Phased Array System, Radar or Communications; `quantiz` is a Communications Toolbox function, so the quantizer is written out). Phased Array System Toolbox functions (`ambgfun`, `phased.LinearFMWaveform`) are not needed; if a licence becomes available they are a cross-check, not a dependency.

## D8 (2026-11-22): mid-rise quantizer by default

The default quantizer is mid-rise (no zero output level, so the sign survives even when the signal is far smaller than one step). A mid-tread quantizer (with a zero level) is run as a comparison in A2, because the two behave very differently when a strong echo has set the gain.

## D9 (2026-11-22): randomised threshold for statistics with few values

Found in the first bit-depth run (E-A1a in `EXPERIMENTS.md`): with one bit and the binary Barker 13 code, the detection statistic can take only a few dozen values, so the 10⁻⁴ quantile falls on a value shared by many noise samples. A plain "greater than" threshold can then give a false-alarm rate below the design value, which would make the receiver look worse than it is. The 13.8 dB loss of Barker 13 at one bit in the first run (`reports/adc/superseded/a1_plain_threshold_losses.csv`) was first suspected to come from this; the repeated run with the randomised threshold gave the same 13.8 dB, so that loss is real (the echo-phase effect, E-PH), and the change only guarantees that every receiver is compared at the same false-alarm rate. The fix is the textbook Neyman-Pearson test for discrete statistics: decide "target" above the threshold, and with probability γ when the statistic equals it, with γ set from the noise-only statistics so the false-alarm rate is exactly the design value (`detection.Threshold`, test `test_randomised_threshold_hits_pfa_for_one_bit_barker`). For continuous statistics γ changes nothing measurable. Every A run was repeated with it. Rejected: reporting each receiver at its own achieved Pfa (receivers would no longer be compared at the same false-alarm rate).

## D10 (2026-11-22): Doppler tolerance read at several Doppler values, not only the first crossing

P4 and Frank lose 3.9 dB at νT = 0.5 (up to νT ≈ 0.5 no other delay matches better, so the peak follows the zero-delay cut, |sinc(0.5)| = 2/π for a constant envelope) but recover to within 0.14 dB (P4 64) at νT = 1, one range cell further on, as an LFM does (`reports/ambiguity/doppler_cuts.csv`). The first 1 dB and 3 dB crossings alone would rank them with Barker 13, which never recovers. The table therefore also gives the peak loss at νT = 0.5, 1, 2 and 4 and the peak's range shift at νT = 1 and 4.

## D11 (2026-11-22): mismatched filters designed and scored at one sample per chip

Least-squares and minimum-PSL filters for phase codes are designed and measured at one sample per chip, the form in which the literature reports them, so the numbers can be compared with published ones. LFM windows are measured at 8 samples per cell, because an LFM has no chips and its sidelobe peaks fall between cells.

## D12 (2026-11-22): report "never reached" instead of a number

A receiver whose Pd never reaches 0.9 by 60 dB is reported as "not reached" (inf in the CSV), not dropped. Found: one bit with the dither removed after the ADC at BT 13 levels off below Pd 0.9 (at 60 dB: 0.85 for LFM, 0.86 for P4, 0.66 for OFDM; `reports/adc/a1_pd_curves.csv`). The subtracted dither adds a noise term that does not shrink as the echo grows, while the one-bit echo itself has a fixed size, so short codes cannot climb above it.

## D13 (2026-11-22): state the run-to-run uncertainty from a replicate, not only from the bootstrap

The paired bootstrap resamples detection trials but keeps each receiver's threshold fixed, so it leaves out the error of the threshold estimated from 2 × 10⁶ noise statistics. A full replicate of the bit-depth run (new seed) moved losses by up to 0.08 dB, outside the combined intervals for 2 of 13 receivers (`reports/adc/a1_replicate_check.csv`). Results are therefore quoted with their bootstrap intervals and the note that a rerun may move them by up to about 0.08 dB. Rejected for now: bootstrapping the threshold too (resampling the noise records for every replicate), which would need the detection statistics stored instead of the decisions; listed as an improvement.
