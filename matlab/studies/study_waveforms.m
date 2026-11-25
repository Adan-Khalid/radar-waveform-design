%STUDY_WAVEFORMS Pulse compression metrics of every waveform at BT 13 and BT 64, and a
%   figure of the compressed pulses. Writes reports/matlab/waveform_metrics.csv and
%   reports/matlab/compressed_pulses_matlab.png.

root = radarwave.repoRoot();
outDir = fullfile(root, 'reports', 'matlab');
if ~isfolder(outDir), mkdir(outDir); end
spc = 8;                                  % samples per range cell (1/B)
designs = {13, {'lfm', 'barker', 'p4', 'ofdm'}; 64, {'lfm', 'frank', 'p4', 'ofdm'}};

rows = {};
fig = figure('Position', [100 100 1200 520], 'Color', 'w');
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact');
for d = 1:2
    bt = designs{d, 1};
    for w = 1:4
        name = designs{d, 2}{w};
        s = radarwave.makeWaveform(name, bt, spc);
        m = radarwave.pulseMetrics(s, [], spc);
        label = radarwave.waveformLabel(name, bt);
        rows(end+1, :) = {bt, label, m.psl_db, m.isl_db, m.width_3db_cells, radarwave.paprDb(s)}; %#ok<SAGROW>
        ax = nexttile(tl);
        [y, lag] = radarwave.correlate(s, s);
        plot(ax, lag / spc, max(20 * log10(abs(y)), -80), 'LineWidth', 1);
        xlim(ax, [-bt bt]); ylim(ax, [-60 3]); grid(ax, 'on');
        title(ax, sprintf('%s: PSL %.1f dB', label, m.psl_db));
        if w == 1, ylabel(ax, sprintf('BT = %d, output (dB)', bt)); end
        if d == 2, xlabel(ax, 'delay (range cells)'); end
    end
end
title(tl, 'Compressed pulses at equal bandwidth and duration (MATLAB)');
T = cell2table(rows, 'VariableNames', {'time_bandwidth', 'waveform', 'psl_db', 'isl_db', ...
    'width_3db_cells', 'papr_db'});
writetable(T, fullfile(outDir, 'waveform_metrics.csv'));
exportgraphics(fig, fullfile(outDir, 'compressed_pulses_matlab.png'), 'Resolution', 150);
disp(T);
