function stats = noiseStatistics(s, receivers, nStats, seed)
%NOISESTATISTICS Noise-only detection statistics |y|^2 of every receiver (same noise).
%   Long noise records are quantized once and correlated at every lag where the
%   waveform fits, so each record gives many statistics. Returns a cell array.
rng(seed, 'twister');
s = s(:).';
n = numel(s);
m = 16 * n;
perRow = m - n + 1;
rowsTotal = ceil(nStats / perRow);
rowsPerBatch = max(1, min(rowsTotal, floor(2e6 / m)));
sc = 1 / sqrt(2);
stats = cell(numel(receivers), 1);
for k = 1:numel(receivers)
    stats{k} = zeros(rowsTotal * perRow, 1);
end
done = 0;
while done < rowsTotal
    rows = min(rowsPerBatch, rowsTotal - done);
    re = randn(rows, m) * sc;
    im = randn(rows, m) * sc;
    uRe = rand(rows, m) - 0.5;
    uIm = rand(rows, m) - 0.5;
    first = done * perRow + 1;
    last = (done + rows) * perRow;
    for k = 1:numel(receivers)
        z = radarwave.applyReceiver(receivers(k), re, im, uRe, uIm);
        y = radarwave.correlate(z, s, 'valid');
        v = abs(y.').^2;              % (perRow x rows)
        stats{k}(first:last) = v(:);
    end
    done = done + rows;
end
for k = 1:numel(receivers)
    stats{k} = stats{k}(1:nStats);
end
end
