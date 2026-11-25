function cut = dopplerCutSummary(s, ref, nuT, samplesPerCell)
%DOPPLERCUTSUMMARY Peak loss (dB), peak delay (cells) and PSL at each Doppler value.
if nargin < 2, ref = []; end
if nargin < 4, samplesPerCell = 1; end
nuT = nuT(:).';
peak = zeros(size(nuT));
shift = zeros(size(nuT));
psl = zeros(size(nuT));
chunk = 256;
for start = 1:chunk:numel(nuT)
    sel = start:min(start + chunk - 1, numel(nuT));
    [chi, lags] = radarwave.ambiguity(s, ref, nuT(sel));
    mag = abs(chi);
    [pk, k] = max(mag, [], 1);
    peak(sel) = pk;
    shift(sel) = lags(k).' / samplesPerCell;
    for j = 1:numel(sel)
        psl(sel(j)) = radarwave.pslDb(mag(:, j));
    end
end
chi0 = radarwave.ambiguity(s, ref, 0);
cut = struct('nu_t', nuT, 'peak_loss_db', 20 * log10(peak / max(abs(chi0))), ...
    'shift_cells', shift, 'psl_db', psl);
end
