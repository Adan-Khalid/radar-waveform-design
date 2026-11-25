# Research log

Sources checked on 2026-11-22. "Not verified" marks a detail that could not be confirmed on an official or indexed page; such details are not used in the README until confirmed.

## 1. Textbooks and classic papers the code follows

| Topic | Reference | What is used from it |
| --- | --- | --- |
| Waveforms, ambiguity function | N. Levanon and E. Mozeson, *Radar Signals*, Wiley-IEEE Press, 2004 (ISBN 978-0-471-47378-7) | matched filter, ambiguity function, Barker, Frank and P4 codes, multicarrier phase-coded pulses |
| Radar signal processing | M. A. Richards, *Fundamentals of Radar Signal Processing*, 2nd ed., McGraw-Hill, 2014 (ISBN 978-0-07-179832-7) | sampling and quantization, matched filter, windows, detection |
| Ambiguity function origin | P. M. Woodward, *Probability and Information Theory, with Applications to Radar*, Pergamon, 1953, chapter 7 | definition of the ambiguity function |
| Frank code | R. L. Frank, "Polyphase codes with good nonperiodic correlation properties", IEEE Trans. Information Theory, IT-9(1), 43-45, 1963, doi:10.1109/TIT.1963.1057798 | Frank phases |
| P4 code | B. L. Lewis and F. F. Kretschmer, "Linear frequency modulation derived polyphase pulse compression codes", IEEE Trans. Aerospace and Electronic Systems, no. 5, 637-641, 1982, doi:10.1109/TAES.1982.309276 (volume number not shown on the record; month not verified) | P4 phases and their Doppler tolerance |
| Mismatched filter (least squares) | M. H. Ackroyd and F. Ghani, "Optimum mismatched filters for sidelobe suppression", IEEE TAES, AES-9(2), 214-218, 1973, doi:10.1109/TAES.1973.309769 | least-squares receive filter |
| Mismatched filter (peak sidelobe) | J. M. Baden and M. N. Cohen, "Optimal peak sidelobe filters for biphase pulse compression", IEEE International Radar Conference, 249-252, 1990, doi:10.1109/RADAR.1990.201171 | optional minimum-PSL filter |
| Codes and receive filters | P. Stoica, J. Li and M. Xue, "Transmit codes and receive filters for radar", IEEE Signal Processing Magazine, 25(6), 94-109, 2008, doi:10.1109/MSP.2008.929231 | sidelobe level against SNR loss trade-off |
| OFDM radar | C. Sturm and W. Wiesbeck, "Waveform design and signal processing aspects for fusion of wireless communications and radar sensing", Proc. IEEE, 99(7), 1236-1259, 2011, doi:10.1109/JPROC.2011.2131110 | OFDM pulse; the reciprocal filter is left for later work |
| Windows | F. J. Harris, "On the use of windows for harmonic analysis with the discrete Fourier transform", Proc. IEEE, 66(1), 51-83, 1978, doi:10.1109/PROC.1978.10837 | processing loss = 10 log10(ENBW): Hamming 1.34 dB, Hann 1.76 dB (Table I lists ENBW 1.36 and 1.50 bins) |
| Bussgang | J. J. Bussgang, "Crosscorrelation functions of amplitude-distorted Gaussian signals", MIT RLE Technical Report 216, 1952, hdl:1721.1/4847 | a Gaussian input and its quantized version are correlated as gain times input: the basis of the quantizer efficiency |
| Clipped noise | J. H. Van Vleck and D. Middleton, "The spectrum of clipped noise", Proc. IEEE, 54(1), 2-19, 1966, doi:10.1109/PROC.1966.4567 | background on hard limiting |
| Optimal quantizers | J. Max, "Quantizing for minimum distortion", IRE Trans. Information Theory, IT-6(1), 7-12, 1960, doi:10.1109/TIT.1960.1057548 | Table II, best uniform step for a unit Gaussian: 1.596 (2 levels), 0.9957 (4), 0.5860 (8), 0.3352 (16) |
| Dither | R. M. Gray and T. G. Stockham, "Dithered quantizers", IEEE Trans. Information Theory, 39(3), 805-812, 1993, doi:10.1109/18.256489; S. P. Lipshitz, R. A. Wannamaker and J. Vanderkooy, "Quantization and dither: a theoretical survey", J. Audio Eng. Soc., 40(5), 355-375, 1992 | subtractive and non-subtractive dither properties |
| Hard limiter SNR | W. B. Davenport, "Signal-to-noise ratios in band-pass limiters", J. Applied Physics, 24(6), 720-727, 1953, doi:10.1063/1.1721365 | one signal in noise: output/input SNR tends to π/4 at low SNR. Two signals: J. J. Jones, "Hard-limiting of two signals in random noise", IEEE Trans. IT, 9(1), 34-42, 1963, doi:10.1109/TIT.1963.1057797. The often quoted "6 dB weak-signal suppression" was not verified and is not used |
| Required SNR | W. J. Albersheim, Proc. IEEE, 69(7), 839, 1981, doi:10.1109/PROC.1981.12082; D. A. Shnidman, "Determination of required SNR values", IEEE TAES, 38(3), 1059-1064, 2002, doi:10.1109/TAES.2002.1039422 | approximations to compare with the exact Marcum Q value used here |

## 2. Low-resolution ADCs in radar and sensing: closest work

| Reference | What it does | Limit or difference from this project |
| --- | --- | --- |
| M. A. Shaikh, G. Joseph, A. Pandharipande and N. J. Myers, "Situation-aware dithering for dynamic range enhancement of a mixed-resolution ADC in automotive radar receivers", IEEE Trans. Radar Systems, 3, 1474-1488, 2025, doi:10.1109/TRS.2025.3620087 (conference version: IEEE RadarConf 2025, doi:10.1109/RadarConf2559087.2025.11205063) | a high-resolution first slot steers the dither of later one-bit slots, so strong targets stop masking weak ones | FMCW and a mixed-resolution design; this project is a controlled single-pulse comparison of waveforms, bit depths and dither types against theory |
| A. Ameri, A. Bose, J. Li and M. Soltanalian, "One-bit radar processing with time-varying sampling thresholds", IEEE Trans. Signal Processing, 67(20), 5297-5308, 2019, doi:10.1109/TSP.2019.2939086 | one-bit radar with time-varying thresholds (a form of dither), recovered by weighted least squares | estimation algorithms, not detection loss against bit depth |
| B. Jin, J. Zhu, Q. Wu, Y. Zhang and Z. Xu, "One-bit LFMCW radar: spectrum analysis and target detection", IEEE TAES, 56(4), 2732-2750, 2020 (arXiv 1905.09440; journal DOI not verified) | one-bit harmonics cause false alarms; about 1 dB loss reported with a strong target | FMCW, one bit only |
| T. Feuillen, C. Xu, J. Louveaux, L. Vandendorpe and L. Jacques, "Quantity over quality: dithered quantization for compressive radar systems", IEEE RadarConf 2019 (arXiv 1811.11575; DOI not verified) | dither removes range ambiguities down to one bit | compressive radar setting |
| F. Xi, N. Shlezinger and Y. C. Eldar, "BiLiMO: bit-limited MIMO radar via task-based quantization", IEEE Trans. Signal Processing, 2021, doi:10.1109/TSP.2021.3116343 | analog combining before few-bit ADCs | MIMO, task-based design |
| C. Studer and G. Durisi, "Quantized massive MU-MIMO-OFDM uplink", IEEE Trans. Communications, 64(6), 2387-2399, 2016, doi:10.1109/TCOMM.2016.2558151 | about 4 to 6 bits cost little in a communication link | communication, not radar detection |
| H. A. Harms, J. Palmer, S. J. Searle and L. M. Davis, "Impact of quantization on passive radar target detection", IET Radar 2012 | ADC bits and dynamic range for OFDM passive radar with a strong direct signal | numeric results not verified |
| H. Xiao, S. Yang and W. Yi, arXiv 2106.00612, 2021 | 3-bit Rao-test detection close to unquantized | preprint |
| J. Zhu, H. Zhang, N. Zhang, J. Fang and F. Qu, arXiv 2307.00491 (v4, 2024) | detection of weak sinusoids with few-bit ADCs next to strong ones | line spectra, not pulse compression |

## 3. Existing tools

| Tool | What it offers | Why it is not used as the core here |
| --- | --- | --- |
| MATLAB `ambgfun`, `phased.LinearFMWaveform`, `phased.PhaseCodedWaveform` (Phased Array System Toolbox) | ambiguity functions; LFM with rectangular or Gaussian envelope; Frank, P1 to P4, Zadoff-Chu, Barker (up to 13) | not in the free MATLAB Online licence; usable later as a cross-check |
| `sdr` (Python, mhostetter) | Barker and Zadoff-Chu codes, replica correlator, Albersheim and Shnidman | no ambiguity function |
| RadarSimPy (GPL-3.0) | radar scene simulator with a C++ engine | no ambiguity function, pulse compression metrics or ADC bit study in its README |
| HermesPy (AGPL) | sensing and communication link simulator | no ambiguity function or quantization study in its README |
| GRAF (research code, 2025) | differentiable ambiguity function for waveform optimisation | research prototype |

No maintained Python package for the narrowband ambiguity function was found, so it is written here (and in MATLAB), with tests.

## 4. MATLAB Online (free licence), checked on the MathWorks page

- Products: MATLAB, Simulink, Control System, Curve Fitting, Deep Learning, Image Processing, Optimization, Signal Processing, Statistics and Machine Learning, Symbolic Math and Text Analytics toolboxes. No Phased Array System, Radar or Communications Toolbox.
- 20 hours per month, 15 minutes of continuous compute, 15 minutes idle timeout, 5 GB MATLAB Drive.
- Functions: `taylorwin`, `chebwin`, `hamming`, `hann`, `marcumq` are Signal Processing Toolbox; `xcorr` and `besseli` are base MATLAB; `ncx2cdf` and `binofit` are Statistics; `quantiz` is Communications Toolbox (not available, so the quantizer is written out); `linprog` is Optimization Toolbox (available, so the optional minimum-PSL filter can run in MATLAB too).
- Consequence: MATLAB Monte Carlo runs use fewer trials and are split so each call finishes well under 15 minutes.

## 5. Theory pre-checks (scratch calculations before writing the package)

Kept as a record of what was expected before the code existed. The final values are the ones in `reports/` (for example the OFDM 64 median PAPR is 7.6 dB at 8 samples per cell and 6.8 dB at one sample per cell in `reports/waveforms/ofdm_spread.csv`).

These were computed with short scratch scripts to set test tolerances and trial counts. The repository's own scripts will recompute every one of them into `reports/`; until then they are not results.

- Quantizer efficiency for a Gaussian input at the best uniform step (Bussgang): loss 1.961, 0.550, 0.166, 0.050 dB for 1, 2, 3, 4 bits, at steps 0.9958, 0.5860, 0.3352 σ (2, 3, 4 bits), matching Max's Table II. At those steps, non-subtractive dither of one step peak to peak raises the loss to 0.90, 0.29 and 0.09 dB: when noise dominates, dither only adds noise.
- Required SNR for Pd 0.9, single pulse, non-fluctuating target, Pfa 10⁻⁴: 11.75 dB (Marcum Q); 10.76 dB at Pfa 10⁻³.
- A quick Monte Carlo (4,000 trials, Pfa 10⁻³, noisy) gave a 1-bit loss of 3.9, 2.3, 2.0 and 2.0 dB at BT 13, 64, 128 and 1024: the 1.96 dB theory needs a low SNR per sample, so short codes lose more. Run A3 measures this properly.
- Sidelobes at 8 samples per 1/B: Barker 13 PSL −22.28 dB, ISL −11.49 dB (exact); LFM PSL −14.4, −13.5 and −13.4 dB at BT 13, 64 and 1000 (ISL −10.7, −9.8, −9.7 dB); P4 13 −17.3 dB; P4 64 −24.4 dB; Frank 64 −27.8 dB; OFDM with random QPSK median PSL −12.6 dB (13 subcarriers) and −13.5 dB (64), median PAPR 5.7 and 7.3 dB.
- Window SNR losses: Hann 1.77 dB, Hamming 1.35 dB, Blackman 2.38 dB, Taylor (35 dB, n̄ = 5) 0.93 dB.
- Speed: with float32 and the same noise shared by all receivers, about 7 x 10⁵ detection statistics per second at BT 64 on this laptop, so the full bit-depth run takes minutes.

## 6. Baselines

- Sidelobe control: the plain matched filter (no window, no mismatch) is the baseline for every window and mismatched filter.
- ADC study: the unquantized receiver is the baseline; the Bussgang efficiency and the Marcum Q curve are the theory references.

## 7. Candidate research question

For pulse waveforms of equal time-bandwidth product, how does the detection loss of a 1 to 4-bit receiver depend on how the ADC full scale is set (by the noise, or by a strong echo) and on the waveform's peak-to-average power ratio, and when does dither (non-subtractive or subtractive) reduce the loss rather than add noise?

Status: candidate only. Each piece is known (the 2/π one-bit loss, dither theory, strong targets masking weak ones in one-bit radar and dithering as the fix in Shaikh et al. 2025). What is not found in one place is a controlled, reproducible comparison across LFM, phase-coded and OFDM pulses, measured against the Bussgang prediction with confidence intervals. The expectation from theory is that the waveform does not matter while noise sets the full scale, and starts to matter once a strong echo or the OFDM peak-to-average ratio does; the runs test that.

## 8. Open items

- Journal DOIs not verified: Jin et al. 2020, Feuillen et al. 2019, Myers and Heath 2019 (IEEE TWC, doi:10.1109/TWC.2019.2909865; volume and pages not verified).
- Levanon and Mozeson: book DOI not verified (ISBN used instead).
