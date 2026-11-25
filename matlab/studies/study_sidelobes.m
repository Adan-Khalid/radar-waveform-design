%STUDY_SIDELOBES Windows on LFM and mismatched filters on phase codes, with SNR cost.
%   Writes reports/matlab/lfm_windows.csv and mismatched_filters.csv.

root = radarwave.repoRoot();
outDir = fullfile(root, 'reports', 'matlab');
if ~isfolder(outDir), mkdir(outDir); end
spc = 8;

windows = {'rectangular', 'hann', 'hamming', 'taylor', 'blackman'};
rows = {};
for bt = [13 64]
    s = radarwave.lfm(bt, spc);
    for k = 1:numel(windows)
        w = radarwave.windowFor(windows{k}, numel(s), 35, 5);
        m = radarwave.pulseMetrics(s, s .* w, spc);
        rows(end+1, :) = {bt, string(windows{k}), m.psl_db, m.isl_db, m.width_3db_cells, ...
            m.mismatch_loss_db, radarwave.windowLossDb(w)}; %#ok<SAGROW>
    end
end
Tw = cell2table(rows, 'VariableNames', {'time_bandwidth', 'window', 'psl_db', 'isl_db', ...
    'width_3db_cells', 'mismatch_loss_db', 'window_loss_formula_db'});
writetable(Tw, fullfile(outDir, 'lfm_windows.csv'));
disp(Tw);

codes = {'barker', 13; 'p4', 13; 'frank', 64; 'p4', 64};
rows = {};
for c = 1:size(codes, 1)
    n = codes{c, 2};
    code = radarwave.makeWaveform(codes{c, 1}, n, 1);
    label = radarwave.waveformLabel(codes{c, 1}, n);
    m = radarwave.pulseMetrics(code, code, 1);
    rows(end+1, :) = {label, "matched", 1, m.psl_db, m.isl_db, m.mismatch_loss_db}; %#ok<SAGROW>
    for f = 1:4
        ref = radarwave.lsMismatchedFilter(code, f * n, 1e-9);
        m = radarwave.pulseMetrics(code, ref, 1);
        rows(end+1, :) = {label, "least squares", f, m.psl_db, m.isl_db, m.mismatch_loss_db}; %#ok<SAGROW>
    end
    if strcmp(codes{c, 1}, 'barker')
        ref = radarwave.minPslFilterReal(code, 3 * n);
        m = radarwave.pulseMetrics(code, ref, 1);
        rows(end+1, :) = {label, "minimum PSL (LP)", 3, m.psl_db, m.isl_db, m.mismatch_loss_db}; %#ok<SAGROW>
    end
end
Tm = cell2table(rows, 'VariableNames', {'code', 'filter', 'length_factor', 'psl_db', 'isl_db', ...
    'mismatch_loss_db'});
writetable(Tm, fullfile(outDir, 'mismatched_filters.csv'));
disp(Tm);
