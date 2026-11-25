%STUDY_CROSSCHECK Compare MATLAB with the Python results in reports/.
%   Deterministic numbers must agree to a small tolerance (the CSV files keep 8
%   significant digits); Monte Carlo losses must agree within their 95 % intervals.
%   Prints PASS or FAIL per check and writes reports/matlab/crosscheck.csv.

root = radarwave.repoRoot();
rep = fullfile(root, 'reports');
outDir = fullfile(rep, 'matlab');
if ~isfolder(outDir), mkdir(outDir); end
results = {};
check = @(name, diff, tol) {name, diff, tol, diff <= tol};

% 1. Waveform samples against Python's (exported with 17 digits)
S = readtable(fullfile(rep, 'crosscheck', 'waveform_samples.csv'), 'TextType', 'string');
cases = {"lfm_bt13_spc8", radarwave.lfm(13, 8); ...
         "p4_13_spc8", radarwave.phaseCode(radarwave.p4Phases(13), 8); ...
         "frank_64_spc1", radarwave.phaseCode(radarwave.frankPhases(64), 1); ...
         "ofdm_13_spc2", radarwave.ofdm(radarwave.ofdmSymbols(13), 2)};
for c = 1:size(cases, 1)
    sel = S.signal == cases{c, 1};
    py = complex(S.re(sel), S.im(sel));
    results(end+1, :) = check("samples " + cases{c, 1}, max(abs(py - cases{c, 2})), 1e-12); %#ok<SAGROW>
end

% 2. Pulse compression metrics, every waveform, both design points
P = readtable(fullfile(rep, 'waveforms', 'metrics.csv'), 'TextType', 'string');
designs = {13, {'lfm', 'barker', 'p4', 'ofdm'}; 64, {'lfm', 'frank', 'p4', 'ofdm'}};
worst = 0;
for d = 1:2
    bt = designs{d, 1};
    for w = 1:4
        name = designs{d, 2}{w};
        s = radarwave.makeWaveform(name, bt, 8);
        m = radarwave.pulseMetrics(s, [], 8);
        row = P(P.waveform == radarwave.waveformLabel(name, bt), :);
        worst = max([worst, abs(m.psl_db - row.psl_db), abs(m.isl_db - row.isl_db), ...
            abs(m.width_3db_cells - row.width_3db_cells), abs(radarwave.paprDb(s) - row.papr_db)]);
    end
end
results(end+1, :) = check("waveform metrics (PSL, ISL, width, PAPR)", worst, 1e-5);

% 3. Least-squares filter taps, Barker 13, 39 taps
F = readtable(fullfile(rep, 'crosscheck', 'ls_filter_barker13_x3.csv'));
ref = radarwave.lsMismatchedFilter(radarwave.phaseCode(radarwave.barkerPhases(13)), 39, 1e-9);
results(end+1, :) = check("LS filter taps, Barker 13 x3", max(abs(complex(F.re, F.im) - ref)), 1e-9);

% 4. Windows on LFM and mismatched filters (same row order as Python)
W = readtable(fullfile(rep, 'sidelobes', 'lfm_windows.csv'), 'TextType', 'string');
windows = {'rectangular', 'hann', 'hamming', 'taylor', 'blackman'};
worst = 0; r = 0;
for bt = [13 64]
    s = radarwave.lfm(bt, 8);
    for k = 1:numel(windows)
        r = r + 1;
        w = radarwave.windowFor(windows{k}, numel(s), 35, 5);
        m = radarwave.pulseMetrics(s, s .* w, 8);
        worst = max([worst, abs(m.psl_db - W.psl_db(r)), abs(m.isl_db - W.isl_db(r)), ...
            abs(m.mismatch_loss_db - W.mismatch_loss_db(r))]);
    end
end
results(end+1, :) = check("LFM windows (PSL, ISL, loss)", worst, 1e-5);

M = readtable(fullfile(rep, 'sidelobes', 'mismatched_filters.csv'), 'TextType', 'string');
codes = {'barker', 13; 'p4', 13; 'frank', 64; 'p4', 64};
worst = 0; worstLp = 0;
for c = 1:size(codes, 1)
    n = codes{c, 2};
    code = radarwave.makeWaveform(codes{c, 1}, n, 1);
    label = radarwave.waveformLabel(codes{c, 1}, n);
    for f = 1:4
        m = radarwave.pulseMetrics(code, radarwave.lsMismatchedFilter(code, f * n, 1e-9), 1);
        row = M(M.code == label & M.filter == "least squares" & M.length_factor == f, :);
        worst = max([worst, abs(m.psl_db - row.psl_db), abs(m.mismatch_loss_db - row.mismatch_loss_db)]);
    end
    if strcmp(codes{c, 1}, 'barker')
        m = radarwave.pulseMetrics(code, radarwave.minPslFilterReal(code, 3 * n), 1);
        row = M(M.code == label & startsWith(M.filter, "minimum"), :);
        worstLp = abs(m.psl_db - row.psl_db);   % LP optimum value; solvers may pick other taps
    end
end
results(end+1, :) = check("least-squares filters (PSL, loss)", worst, 1e-5);
results(end+1, :) = check("minimum-PSL filter, Barker 13 (PSL)", worstLp, 0.05);

% 5. Doppler tolerance (needs study_ambiguity first)
fileM = fullfile(outDir, 'doppler_tolerance.csv');
if isfile(fileM)
    A = readtable(fileM, 'TextType', 'string');
    B = readtable(fullfile(rep, 'ambiguity', 'doppler_tolerance.csv'), 'TextType', 'string');
    worst = 0;
    for k = 1:height(A)
        row = B(B.waveform == A.waveform(k), :);
        worst = max([worst, abs(A.nu_t_1db(k) - row.nu_t_1db), abs(A.nu_t_3db(k) - row.nu_t_3db), ...
            abs(A.loss_at_1_db(k) - row.loss_at_1_db), abs(A.loss_at_4_db(k) - row.loss_at_4_db)]);
    end
    results(end+1, :) = check("Doppler tolerance", worst, 1e-5);
end

% 6. Theory: required SNR and quantizer losses
R = readtable(fullfile(rep, 'theory', 'required_snr.csv'));
worst = 0;
for k = 1:height(R)
    worst = max(worst, abs(radarwave.requiredSnrDb(R.pd(k), R.pfa(k)) - R.required_snr_marcum_db(k)));
end
results(end+1, :) = check("required SNR (Marcum Q)", worst, 1e-5);
Q = readtable(fullfile(rep, 'theory', 'quantizer_loss_optimal.csv'), 'TextType', 'string');
worst = 0;
for k = 1:height(Q)
    eta = radarwave.quantizerEfficiency(Q.bits(k), Q.a_over_sigma(k), char(Q.kind(k)), char(Q.dither(k)));
    worst = max(worst, abs(-10 * log10(eta) - Q.loss_db(k)));
end
results(end+1, :) = check("quantizer losses (Bussgang)", worst, 1e-5);

% 7. Monte Carlo: MATLAB ADC losses inside the Python and MATLAB intervals combined
fileA = fullfile(outDir, 'adc_losses.csv');
fileP = fullfile(rep, 'adc', 'a1_losses.csv');
if isfile(fileA) && isfile(fileP)
    A = readtable(fileA, 'TextType', 'string');
    B = readtable(fileP, 'TextType', 'string');
    B = B(B.waveform == "LFM (BT 64)", :);
    for k = 2:height(A)
        row = B(B.receiver == A.receiver(k), :);
        if isempty(row), continue; end
        halfM = (A.loss_hi(k) - A.loss_lo(k)) / 2;
        halfP = (row.loss_hi - row.loss_lo) / 2;
        results(end+1, :) = check("ADC loss, " + A.receiver(k), abs(A.loss_db(k) - row.loss_db), ...
            sqrt(halfM^2 + halfP^2)); %#ok<SAGROW>
    end
end

T = cell2table(results, 'VariableNames', {'check', 'difference', 'tolerance', 'pass'});
writetable(T, fullfile(outDir, 'crosscheck.csv'));
for k = 1:height(T)
    if T.pass(k), verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('%s  %-45s difference %.3g (tolerance %.3g)\n', verdict, T.check(k), ...
        T.difference(k), T.tolerance(k));
end
