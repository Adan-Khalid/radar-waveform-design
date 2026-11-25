%STUDY_AMBIGUITY Ambiguity functions (BT 64) and Doppler tolerance of every waveform.
%   Writes reports/matlab/doppler_tolerance.csv and ambiguity_functions_matlab.png.

root = radarwave.repoRoot();
outDir = fullfile(root, 'reports', 'matlab');
if ~isfolder(outDir), mkdir(outDir); end
spc = 8;
designs = {13, {'lfm', 'barker', 'p4', 'ofdm'}; 64, {'lfm', 'frank', 'p4', 'ofdm'}};

rows = {};
for d = 1:2
    bt = designs{d, 1};
    for w = 1:4
        name = designs{d, 2}{w};
        s = radarwave.makeWaveform(name, bt, spc);
        summ = radarwave.dopplerTolerance(s, [], 4, 0.01, spc);
        if ~isfinite(summ.nu_t_1db) || ~isfinite(summ.nu_t_3db)
            wide = round((0:0.05:bt) * 1e10) / 1e10;   % Doppler-tolerant pulse: scan to BT
            cut2 = radarwave.dopplerCutSummary(s, [], wide, spc);
            if ~isfinite(summ.nu_t_1db)
                summ.nu_t_1db = radarwave.firstCrossing(wide, cut2.peak_loss_db, -1);
            end
            if ~isfinite(summ.nu_t_3db)
                summ.nu_t_3db = radarwave.firstCrossing(wide, cut2.peak_loss_db, -3);
            end
        end
        rows(end+1, :) = {bt, radarwave.waveformLabel(name, bt), summ.nu_t_1db, summ.nu_t_3db, ...
            summ.loss_at_0p5_db, summ.loss_at_1_db, summ.loss_at_2_db, summ.loss_at_4_db, ...
            summ.shift_at_1_cells}; %#ok<SAGROW>
        fprintf('%-12s 3 dB at nuT %.3f, loss at nuT 1: %.2f dB\n', ...
            radarwave.waveformLabel(name, bt), summ.nu_t_3db, summ.loss_at_1_db);
    end
end
T = cell2table(rows, 'VariableNames', {'time_bandwidth', 'waveform', 'nu_t_1db', 'nu_t_3db', ...
    'loss_at_0p5_db', 'loss_at_1_db', 'loss_at_2_db', 'loss_at_4_db', 'shift_at_1_cells'});
writetable(T, fullfile(outDir, 'doppler_tolerance.csv'));

% Ambiguity function images, BT 64, delay +-8 cells, nu*T +-4
fig = figure('Position', [100 100 1200 360], 'Color', 'w');
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact');
spcImg = 4;
nu = linspace(-4, 4, 241);
lags = -8*spcImg:8*spcImg;
for w = 1:4
    name = designs{2, 2}{w};
    s = radarwave.makeWaveform(name, 64, spcImg);
    chi = radarwave.ambiguity(s, [], nu, lags);
    ax = nexttile(tl);
    imagesc(ax, lags / spcImg, nu, max(20 * log10(abs(chi.')), -30));
    axis(ax, 'xy'); clim(ax, [-30 0]); colormap(ax, flipud(gray));
    title(ax, radarwave.waveformLabel(name, 64)); xlabel(ax, 'delay (range cells)');
    if w == 1, ylabel(ax, '\nu T'); end
end
cb = colorbar(ax); cb.Label.String = '|\chi| (dB)';
title(tl, 'Ambiguity functions, BT = 64 (MATLAB)');
exportgraphics(fig, fullfile(outDir, 'ambiguity_functions_matlab.png'), 'Resolution', 150);
disp(T);
