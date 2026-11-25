# Learning notes

Notes on every part of the project: the idea, the maths, the design choice, the alternatives that were rejected, and three interview questions with short answers. Code references point at `src/radarwave/` (Python) and `matlab/+radarwave/` (MATLAB), which are written to implement the same equations (the MATLAB code has not been run yet). Measured numbers live in `reports/`; the values quoted here are the theory they are checked against.

Symbols: s[n] the transmitted pulse (complex baseband samples, unit energy Σ|s|² = 1), B bandwidth, T duration, BT time-bandwidth product, f_s sample rate, σ² complex noise power, σ_c = σ/√2 noise standard deviation of I or Q, a echo amplitude, ν Doppler shift.

---

## 1. The matched filter and its SNR gain

**Code:** `matched.py` (`correlate`, `matched_filter`); test `test_snr_gain_equals_time_bandwidth`.

### The idea

To find a known pulse in noise, correlate the received samples with a copy of the pulse. At the right delay every sample adds in phase and the noise adds randomly, so the pulse's energy is collected into one output sample.

### The maths

Received: r[n] = a e^{jφ} s[n − d] + w[n], w white with power σ². The filter output at lag k is y[k] = Σ_n r[n + k] s*[n].

- At k = d: signal part a e^{jφ} Σ|s|² = a e^{jφ} (unit energy), noise part Σ w s*, whose variance is σ² Σ|s|² = σ². Output SNR = |a|²/σ².
- **Why this is the best linear filter.** For any filter h, by the Cauchy-Schwarz inequality |Σ h* s|² ≤ Σ|h|² Σ|s|², so the output SNR |a|² |Σ h* s|² / (σ² Σ|h|²) is at most |a|²/σ², with equality only for h ∝ s.
- **SNR gain.** Per sample, the echo power is |a|² |s[n]|² = |a|²/N for a constant-envelope pulse of N samples. Output over input SNR = N. Sampled at f_s = B, N = BT, so the gain is 10 log10(BT): 11.1 dB for BT = 13, 18.1 dB for BT = 64.
- **Range resolution.** The compressed pulse has a mainlobe about 1/B wide whatever T is. Long pulses give energy (detection range), wide bandwidth gives resolution; pulse compression decouples the two. In metres: ΔR = c/(2B); 150 m at B = 1 MHz.

### Design choice

The correlation is done by FFT (zero-padded to a power of two, so it is linear, not circular). Output index k means "echo starts k samples into the window", the convention used throughout.

### Alternatives rejected

- Direct time-domain convolution: same numbers, much slower for long records.
- `numpy.correlate` / MATLAB `xcorr`: fine for one pulse, but their lag conventions differ between languages; one own function keeps Python and MATLAB identical.

### Interview questions

1. *Why does a matched filter maximise SNR?* Cauchy-Schwarz: the output SNR is a normalised inner product of the filter with the pulse, largest when they are parallel.
2. *What is the SNR gain of pulse compression and what sets it?* The time-bandwidth product BT (samples at f_s = B): the filter adds N samples coherently while noise adds in power.
3. *Does a longer pulse improve range resolution?* No. Resolution is about c/2B; a longer pulse at the same bandwidth only adds energy (and a longer blind zone).

---

## 2. Waveforms at equal bandwidth and duration

**Code:** `waveforms.py`; tests in `test_waveforms_metrics.py`; settings `configs/waveforms.yaml`.

### The pulses

- **LFM:** s(t) = exp(jπ(B/T)(t − T/2)²); the instantaneous frequency (1/2π) dφ/dt = (B/T)(t − T/2) sweeps −B/2 to +B/2.
- **Barker 13:** 13 chips with phases 0 or π (+ + + + + − − + + − + − +). Its aperiodic autocorrelation has peak 13 and every sidelobe 0 or 1, the best a binary code can do: PSL = 20 log10(1/13) = −22.3 dB.
- **Frank (N = M²):** chip (i, k) has phase 2π i k / M. It is a stepped approximation of LFM (M frequency steps of M chips each).
- **P4 (any N):** φ_n = π n²/N − π n, the phases of an LFM sampled once per chip, centred so the phase steps are symmetric. More Doppler tolerant than Frank because the "frequency" changes smoothly, not in steps.
- **OFDM:** N subcarriers at frequencies (k − (N − 1)/2)/T carrying random QPSK symbols. It is the waveform of a communication link used as a radar pulse (joint sensing and communication).

Each phase code has chip length T/N, so its nominal bandwidth is B = N/T and BT = N; OFDM has N subcarriers 1/T apart, so B = N/T as well. That is how the comparison is made fair: at BT = 13 the four pulses have the same B, the same T, the same energy.

### Design choice

Two design points (BT 13 and 64; `docs/DECISIONS.md` D1). Waveform metrics use 8 samples per 1/B so peaks between range cells are not missed; the ADC study uses 1 sample per cell (D4).

### Alternatives rejected

- Matching a 3 dB or 99 % power bandwidth instead of the nominal one: rectangular-chip codes have sinc spectra, so this gives every code a different chip length and loses the textbook check values.
- Low-PAPR OFDM symbols (Newman phases): a better radar pulse but no longer a data-carrying one; left for future work on OFDM and OTFS sensing.

### Interview questions

1. *Why is P4 more Doppler tolerant than Frank?* Both approximate LFM, but P4's phase steps follow a smooth quadratic, so a Doppler shift mostly moves the peak in delay (like LFM) instead of breaking the code into mismatched groups.
2. *Why are there no binary codes better than Barker beyond length 13?* None exist with all sidelobes ≤ 1 (proved for odd lengths above 13 and searched far for even ones); longer binary codes trade to minimum-peak-sidelobe codes with sidelobes of 2 or more.
3. *What does OFDM cost as a radar pulse?* A high peak-to-average power ratio (the amplifier and the ADC must cover the peaks) and noise-like sidelobes that depend on the data symbols.

---

## 3. Measuring a compressed pulse

**Code:** `metrics.py` (`mainlobe_bounds`, `psl_db`, `isl_db`, `width`, `mismatch_loss_db`); tests `test_barker13_sidelobes_exact`, `test_lfm_large_bt_approaches_sinc`, `test_rectangular_chip_mainlobe_width`.

### Definitions

- **Mainlobe:** from the peak out to the first local minimum on each side.
- **PSL:** largest magnitude outside the mainlobe over the peak, in dB.
- **ISL:** energy outside the mainlobe over energy inside it, in dB.
- **Width:** distance between the −3 dB crossings (linear interpolation), in range cells of 1/B.
- **Mismatch loss:** −10 log10( max_k |y_k|² / (‖s‖² ‖h‖²) ) for a receive reference h. Zero for the matched filter, positive otherwise (Cauchy-Schwarz).

### Check values

- Barker 13: PSL −22.28 dB and ISL 10 log10(12/169) = −11.49 dB, at any oversampling (with rectangular chips the output is the code's autocorrelation joined by straight lines, so the ratio of triangle energies equals the ratio of squared values).
- Rectangular chips: the mainlobe is a triangle 1 − |τ|/t_c, so the −3 dB width is 2(1 − 10^(−3/20)) = 0.584 chip.
- LFM at large BT: the compressed pulse tends to a sinc: PSL −13.26 dB, ISL −9.68 dB, width 0.886/B. At BT = 13 the Fresnel ripples of the spectrum change these (see `reports/waveforms/metrics.csv`).

### Interview questions

1. *PSL or ISL: which matters when?* PSL for one strong target hiding a weak neighbour; ISL for many scatterers (clutter) whose sidelobes add up.
2. *Why define the mainlobe at the first minimum and not the first null?* Windowed or oversampled outputs may never reach exactly zero; the first minimum works for every case and equals the first null when there is one.
3. *What does a 0.886/B width mean in metres?* With B = 1 MHz, c/2B = 150 m per cell, so the −3 dB width is about 133 m.

---

## 4. The ambiguity function and Doppler tolerance

**Code:** `ambiguity.py` (`ambiguity`, `doppler_cut_summary`, `doppler_tolerance`); tests `test_volume_identity_exact`, `test_zero_delay_cut_is_a_sampled_sinc`, `test_lfm_ridge_height_and_slope`, `test_doppler_tolerance_phase_code_vs_lfm`.

### The maths

χ(m, ν) = Σ_n s[n] e^{j2πνn/f_s} r*[n − m]: the output at lag m of the filter matched to r when the echo has Doppler ν (Woodward). Computed for every lag by building s[n] r*[n − m] and multiplying by a matrix of Doppler phases (any Doppler grid, not only FFT bins).

Properties used as tests:
- χ(0, 0) = Σ|s|² = 1.
- Zero-Doppler cut = matched-filter output.
- Zero-delay cut of a constant-envelope pulse: |Σ_n (1/N) e^{j2πνTn/N}| = |sin(πνT) / (N sin(πνT/N))|, a sampled sinc with its first null at νT = 1.
- **Volume:** on a full DFT grid of K ≥ N Doppler values, Parseval gives Σ_k |χ(m, k)|² = K Σ_n |s[n]|² |r[n − m]|², and summing over m gives (1/K) Σ_m Σ_k |χ|² = E_s E_r exactly. Sidelobe energy can be moved around the delay-Doppler plane but not removed (the radar uncertainty principle).
- **LFM ridge:** a Doppler shift on a chirp looks like a time shift. Multiplying by e^{j2πνt} turns exp(jπ(B/T)(t − T/2)²) into the same chirp moved by νT/B, so the peak moves by νT range cells (sign set by the sweep direction) and drops to 1 − |ν|/B.

### Doppler tolerance

From the cuts: the νT at which the peak (over delay) has lost 1 dB and 3 dB, the delay shift of the peak per unit νT (range-Doppler coupling), and the PSL at νT = 0.25 and 0.5. A constant-envelope code with no coupling follows |sinc(νT)|: −3 dB at νT = 0.443. This holds for any code whose peak cannot move to a better delay before νT ≈ 0.5. LFM loses 3 dB only near |ν| = 0.29 B, about 0.66 BT times further out (42 times at BT 64), but gives a range error of νT cells.

### Alternatives rejected

- MATLAB `ambgfun`: needs the Phased Array System Toolbox (not in the free licence) and normalises differently; kept as a later cross-check only.
- FFT over a zero-padded Doppler axis: gives the full plane but a fixed grid; the matrix form allows a fine grid exactly where the tolerance is read.

### Interview questions

1. *What does the ambiguity function tell a radar designer?* For every delay and Doppler mismatch, how much of the matched-filter peak survives and where the sidelobes are: resolution, Doppler tolerance and range-Doppler coupling in one picture.
2. *Why is LFM called Doppler tolerant, and what is the catch?* The peak survives large Doppler shifts, but it moves in range by νT cells (coupling), so a fast target appears at a slightly wrong range.
3. *Can a waveform have a thumbtack ambiguity function with no sidelobes anywhere?* No. The volume is fixed (E²); a thumbtack spreads it as a low floor over the whole plane.

---

## 5. Sidelobe control with amplitude windows (LFM)

**Code:** `sidelobes.py` (`window`, `taylor`, `windowed_reference`); `theory.window_loss_db`; tests `test_windows_match_scipy`, `test_window_losses_harris`; results `reports/sidelobes/lfm_windows.csv`.

### The idea

An LFM's spectrum is nearly rectangular, so its compressed pulse is nearly a sinc with −13 dB sidelobes. Weighting the receive reference with a window tapers that spectrum (time maps to frequency along the chirp) and lowers the sidelobes, at the price of a wider mainlobe and some SNR.

### The maths

The receive reference is h[n] = w[n] s[n]. For a constant-envelope pulse the mismatch loss is

loss = −10 log10( (Σ w)² / (N Σ w²) ) = −10 log10( (mean w)² / mean(w²) ),

which is the window's processing loss: Hamming 1.34 dB, Hann 1.76 dB, Blackman 2.37 dB, Taylor (35 dB, n̄ = 5) about 0.93 dB. (Harris lists ENBW, the equivalent noise bandwidth in bins; the processing loss is 10 log10 ENBW: 1.36 bins gives 1.34 dB for Hamming.) The Taylor window is written out from its formula (A = acosh(10^(SLL/20))/π, σ_p² = n̄²/(A² + (n̄ − ½)²), coefficients F_m) in both languages and checked against SciPy.

### Design choice

Window on receive only (the transmitter stays at full power, constant envelope). Taylor is the usual radar choice: near-constant sidelobes at the design level with the smallest loss for that level.

### Alternatives rejected

- Weighting the transmitted pulse: lowers transmitted energy and needs a linear amplifier.
- Non-linear FM (shaping the sweep so the spectrum itself is tapered): no receive loss, but its own design problem; noted as a next step.
- Frequency-domain weighting: equivalent for large BT (stationary phase); at BT = 13 time weighting is the simpler and standard form.

### Interview questions

1. *Why does a window lower sidelobes but widen the mainlobe?* Tapering the spectrum smooths its edges (fewer ripples in the transform) but reduces the effective bandwidth.
2. *Which window would you pick and why?* Taylor: it reaches a chosen sidelobe level with the least SNR loss and mainlobe growth; Hamming is a fixed-level classic.
3. *Why does windowing work worse at small BT?* The time-to-frequency mapping of the chirp is only approximate when BT is small (Fresnel ripples), so the spectrum is not tapered as cleanly.

---

## 6. Sidelobe control with a mismatched filter (phase codes)

**Code:** `sidelobes.py` (`ls_mismatched_filter`, `minpsl_filter_real`, `convolution_matrix`); tests `test_ls_filter_constraint_and_gain`, `test_ls_filter_beats_matched_isl_for_p4`, `test_minpsl_filter_has_lower_peak_sidelobe_than_ls`; results `reports/sidelobes/mismatched_filters.csv`.

### The maths

Write the filter output as z = S g, with S the (N + K − 1) × K convolution matrix of the code and g the K filter taps. Fix the peak, z_c = a g = 1 (a = row c of S), and minimise the energy of every other output, ‖S_rest g‖². With a Lagrange multiplier:

g = R⁻¹ aᴴ / (a R⁻¹ aᴴ),   R = S_restᴴ S_rest + δ I

(Ackroyd and Ghani; δ is a tiny diagonal load for numerical safety). Because the matched filter is the unique SNR-optimal filter, any other g costs SNR; the mismatch loss from section 3 measures it. Longer filters (K = 2N to 4N) push sidelobes lower with only a small extra loss.

The peak-sidelobe version (Baden and Cohen) minimises the largest sidelobe instead: minimise t subject to |z_k| ≤ t for k ≠ c and z_c = 1. For a real code this is a linear programme (solved with HiGHS in SciPy, `linprog` in MATLAB's Optimization Toolbox).

### Design choice

Least squares as the main method: closed form, works for complex codes (P4, Frank), and minimises ISL. Designed and evaluated at one sample per chip, the form in which the literature reports these filters.

### Alternatives rejected

- Inverse filter (divide by the code spectrum): fails where the spectrum has nulls and amplifies noise there.
- Minimum-PSL filter for complex codes: needs second-order cone programming, not available in base SciPy or free MATLAB Online.

### Interview questions

1. *Why does a mismatched filter always lose SNR?* The matched filter is the only filter that reaches the Cauchy-Schwarz bound; anything else is a compromise between peak height and sidelobe energy.
2. *What is the catch with mismatched filters under Doppler?* They are designed for zero Doppler; a Doppler shift breaks the exact cancellation, so the sidelobes rise faster than the matched filter's (`reports/sidelobes/mismatched_doppler.csv`).
3. *ISL-optimal or PSL-optimal?* ISL (least squares) for distributed clutter and for a closed-form, fast design; PSL (linear programme) when one strong target must not hide a weak one.

---

## 7. Quantization: uniform quantizers, full scale, mid-rise and mid-tread

**Code:** `adc.py` (`quantize_real`, `quantize_with_dither`, `quantize`, `staircase`, `gaussian_optimal_full_scale`); tests in `test_adc_detection.py`.

### The model

Each of I and Q goes through the same b-bit uniform quantizer with full scale A (the clipping level) and step Δ = 2A/2^b.

- **Mid-rise:** 2^b levels at ±Δ/2, ±3Δ/2, ...; outputs (⌊x/Δ⌋ + ½)Δ, clipped to the outer levels. No zero level: one bit keeps only the sign.
- **Mid-tread:** 2^b − 1 levels at 0, ±Δ, ...; a signal smaller than Δ/2 gives exactly zero.

The full scale is the gain setting (the AGC). Too small: the signal clips. Too large: the signal uses only the central levels. For a Gaussian input there is a best step for each b (Max, 1960): 0.9957 σ_c, 0.5860 σ_c and 0.3352 σ_c for 2, 3 and 4 bits.

### Design choice

Mid-rise by default (D8) and the full scale set from the noise level, the usual setting for a receiver that waits for weak echoes. The full-scale study (A2) then moves A far from that optimum, which is what happens when a strong echo elsewhere in the window sets the gain.

### Alternatives rejected

- Non-uniform (Lloyd-Max) levels: slightly better for Gaussian noise (4 levels: distortion 0.1175 against 0.1188) but not how ADCs are built.
- Communications Toolbox `quantiz`: not in free MATLAB Online, and an own quantizer keeps both languages identical.

### Interview questions

1. *What is the difference between mid-rise and mid-tread, and why does it matter for weak signals?* Mid-tread has a zero level, so a signal smaller than half a step disappears; mid-rise always keeps the sign.
2. *How do you set the gain in front of a low-bit ADC?* So the step suits the noise (about one σ for 2 bits, 0.59 σ for 3 bits) unless strong echoes must not clip; that conflict is the dynamic-range problem.
3. *What does one extra bit buy?* Half the step, so about 6 dB more range between the strongest signal that does not clip and the noise level that is still resolved.

---

## 8. Why one bit costs 2/π: the Bussgang efficiency of a quantizer

**Code:** `theory.py` (`quantizer_efficiency`, `quantizer_loss_db`); tests `test_quantizer_loss_at_max_steps`, `test_one_bit_is_two_over_pi_for_any_full_scale`, `test_efficiency_against_brute_force_integration`; results `reports/theory/quantizer_loss_optimal.csv`.

### The maths

When each sample is mostly noise (low per-sample SNR), the quantizer input x is nearly Gaussian. Bussgang's theorem: the output q(x) can be written as q(x) = g x + e with gain g = E[x q(x)]/σ² and a distortion e uncorrelated with x. A weak signal passes with gain g; the output noise power is E[q²]. So the SNR is multiplied by the efficiency

η = g² σ² / E[q²] = E[x q(x)]² / (σ² E[q(x)²]),

and the detection loss is −10 log10 η, as long as the matched filter's output is close to Gaussian (it sums N samples).

- **One bit:** E[x sign(x)] = σ √(2/π), E[sign²] = 1, so η = 2/π: a loss of 1.96 dB, whatever the full scale.
- **b bits:** with thresholds t_i and levels y_i, E[x q] = σ Σ y_i (φ(t_{i−1}/σ) − φ(t_i/σ)) and E[q²] = Σ y_i² (Φ(t_i/σ) − Φ(t_{i−1}/σ)), exact sums of normal densities and distribution functions. At Max's steps: 0.55, 0.17 and 0.05 dB for 2, 3 and 4 bits. At an MSE-optimal step E[xq] = E[q²], so η = 1 − (distortion/σ²).

### Where the theory stops

It assumes each sample's SNR is low. For a short code the echo is not small against the noise in each sample at the SNR where detection happens (per-sample SNR = SNR − 10 log10 BT, about +0.6 dB at BT = 13 against −6.3 dB at BT = 64), and one bit loses more than 1.96 dB. Run A3 measures how fast the loss approaches the theory as BT grows.

Two effects found in the runs show how one-bit sampling of a short code goes wrong:

- **The echo phase matters.** I and Q are quantized separately. A real (binary) code such as Barker 13 arriving at phase 0 puts all its energy in I; Q then holds only noise signs, and the statistic (sum over I)² + (sum over Q)² is 13² plus the square of a sum of 13 random signs, which crosses the threshold set for Pfa = 10⁻⁴ only when those noise signs happen to line up. A Barker 13 echo at phase 0 is therefore found with Pd 0.09 at most, even at 40 dB, while at 45° it is found with Pd 0.95 at 16 dB (`reports/adc/phase_pd_curves.csv`). A complex code (P4) spreads its echo over I and Q at every phase and is found at all of them. Averaged over random phases this is why Barker 13 needs 13.8 dB more SNR at one bit. With long codes the sums are nearly Gaussian and the effect disappears.
- **Subtracting dither after a one-bit ADC caps Pd.** With the dither removed, the output is (A/2) sign(x + d) − d: the echo part has a fixed size however strong the echo, while −d adds noise that never shrinks. At BT = 13, Pd at 60 dB is still only 0.85 (LFM), 0.86 (P4) and 0.66 (OFDM), so 0.9 is never reached (`reports/adc/a1_pd_curves.csv`).

Both are reasons why one-bit radar receivers use long codes, many pulses, or extra processing (time-varying thresholds, mixed-resolution channels) rather than relying on the 1.96 dB figure.

### Interview questions

1. *Derive the one-bit loss.* For x ~ N(0, σ²): E[x sign x] = 2 ∫₀^∞ x φ dx = σ√(2/π), E[sign²] = 1, so η = 2/π, −1.96 dB.
2. *Why doesn't the waveform matter in this formula?* At low per-sample SNR the quantizer only sees Gaussian noise; the signal is a small perturbation, so the loss depends on the noise statistics and the quantizer, not on the pulse shape.
3. *When would 1-bit sampling lose much more than 2 dB?* When each sample's SNR is not small (short codes, strong echoes), when a strong echo dominates the sign (weak targets suppressed), or when the gain clips or under-loads a multi-bit quantizer.

---

## 9. Dither: non-subtractive and subtractive

**Code:** `adc.py` (`quantize_with_dither`); `theory.quantizer_efficiency(..., dither=...)`; tests `test_nonsubtractive_dither_makes_the_mean_linear`, `test_subtractive_dither_error_is_uniform_and_independent`, `test_dither_efficiency_against_monte_carlo`.

### The maths

Add d ~ U(−Δ/2, Δ/2) before the quantizer.

- **Non-subtractive** (output q(x + d)): the mean output becomes exactly linear, E_d[q(x + d)] = x, while x + d stays inside the range. The quantizer no longer distorts the signal on average, but the output now carries the dither's power too.
- **Subtractive** (output q(x + d) − d, with d known to the receiver): the total error is uniform with variance Δ²/12 and independent of the input (Schuchman's condition; Gray and Stockham), as long as nothing clips.

The efficiency formula of section 8 still applies with the dither averaged in: E[x z]² / (σ² E[z²]), z the output after dither. Expectations over x stay exact; the average over d uses Gauss-Legendre nodes.

### What to expect

When the noise already spans several steps (full scale set for the noise), the noise is its own dither: extra non-subtractive dither only adds power (2 bits: 0.90 dB instead of 0.55 dB).

The full-scale run (A2) shows the less obvious case. When a strong echo elsewhere in the window sets the gain, the step becomes much larger than the noise around a weak target. A mid-rise quantizer then turns into a one-bit receiver for that target and loses at most the one-bit 2.3 dB. Dither of one step peak to peak adds power in proportion to the step squared, so its loss keeps growing (7.9 dB at 30 dB full scale for 3 bits even with the dither removed). Dither earns its place only where the quantization error is a deterministic function of a strong signal in the same samples (spurs), which is the overlapping case left for run A4, or with a mid-tread quantizer, whose zero level would otherwise swallow the weak target.

### Alternatives rejected

- Triangular (two-step) dither: also makes the error power independent of the signal for non-subtractive dither, at twice the dither power; not needed for the questions here.
- Time-varying one-bit thresholds (Ameri et al. 2019): a form of dither aimed at estimating amplitudes from one-bit data; a different goal from detection.

### Interview questions

1. *Why would anyone add noise to a signal on purpose?* To make the quantization error independent of the signal: spurs and dead zones turn into a flat noise floor that averaging and the matched filter can handle.
2. *Subtractive or non-subtractive?* Subtractive removes the dither's power after the ADC but needs the exact dither values (a known pseudo-random source); non-subtractive is simpler and costs dither power.
3. *When does dither hurt detection?* When the noise is already larger than a step: the noise dithers for free and the added dither only raises the floor.

---

## 10. Detection at fixed Pfa: thresholds, Marcum Q, intervals

**Code:** `detection.py` (`Receiver`, `noise_statistics`, `check_pfa`, `TrialSet`, `detections`, `search_required_snr`, `loss_table`, `wilson`, `crossing_db`); `theory.pd_marcum`, `theory.required_snr_db`; tests in `test_adc_detection.py`.

### The maths

The statistic is |y|², y the filter output at the target's cell.

- Noise only, unquantized: |y|²/σ² is exponential, so Pfa = e^{−T} and T = −ln Pfa (9.21 for 10⁻⁴).
- With a non-fluctuating echo of SNR S and random phase: 2|y|²/σ² is non-central χ² with 2 degrees of freedom and non-centrality 2S, so Pd = Q₁(√(2S), √(2T)), Marcum's Q function. Pd = 0.9 at Pfa = 10⁻⁴ needs S = 11.75 dB.

### How the simulation measures loss

1. **Threshold per receiver** from 2 × 10⁶ noise-only statistics of that same receiver (long noise records, quantized once, correlated at every lag), then the false-alarm rate is checked on fresh noise. After a coarse quantizer the statistic may take only a few values (one bit with a 13-chip binary code gives a few dozen), so the threshold is **randomised** on ties, the textbook Neyman-Pearson test for discrete statistics: decide "target" above T, and with probability γ at T, with γ chosen so the false-alarm rate is exactly the design value (`Threshold` in `detection.py`, D9).
2. **Common random numbers:** all receivers see the same noise, phases and dither draws, and the same noise is reused at every SNR. Differences between receivers then have far less Monte Carlo noise than either result alone.
3. **Coarse-to-fine search:** a coarse grid with few trials finds where Pd crosses 0.9; a fine grid (0.1 dB steps, 20,000 trials) measures it. The required SNR is the linear interpolation at Pd = 0.9.
4. **Intervals:** Wilson 95 % intervals for each Pd (better than p ± 1.96√(p(1−p)/n) near 0 and 1); for required SNR and loss, a bootstrap with Poisson(1) weights on the trials, the same weights for every receiver (paired), 1,000 replicates, 2.5 and 97.5 percentiles.

The threshold's own uncertainty is not in the bootstrap; the fresh-noise Pfa check shows its size (`reports/adc/*_pfa_check.csv`), and a full replicate with new random numbers moved the losses by up to 0.08 dB, a little more than the bootstrap intervals in 2 of 13 cases (`reports/adc/a1_replicate_check.csv`, D13).

### Design choice

A clairvoyant detector at the known cell (D5): it measures what the ADC costs, not CFAR or straddle losses, which are separate topics.

### Alternatives rejected

- Theoretical threshold −ln Pfa for every receiver: wrong after quantization (the noise statistic is no longer exactly exponential), so each receiver gets its own empirical threshold.
- Fitting a curve to Pd(SNR): assumes a shape; local interpolation on a fine grid does not.

### Interview questions

1. *Why set Pfa first and then compare Pd?* That is how radars are run (Neyman-Pearson): the false-alarm rate is a requirement, and detectors are compared at equal Pfa.
2. *Why use common random numbers?* Correlated errors cancel in the difference, so a 0.05 dB loss can be measured with a reasonable number of trials.
3. *What is the Wilson interval and why use it?* A binomial interval from inverting the score test; it stays inside [0, 1] and keeps its coverage near Pd = 0.9 and Pfa = 10⁻⁴, where the normal approximation does not.
