function hits = detectionsAt(s, rx, thr, snrDb, trials, seed)
%DETECTIONSAT Detections (numel(snrDb) x trials logical) at the target's cell.
%   The same seed gives the same noise, phases and dither at every SNR and for every
%   receiver (common random numbers), so receivers can be compared in pairs.
rng(seed, 'twister');
s = s(:).';
n = numel(s);
sc = 1 / sqrt(2);
re = randn(trials, n) * sc;
im = randn(trials, n) * sc;
phase = rand(trials, 1) * 2 * pi;
uRe = rand(trials, n) - 0.5;
uIm = rand(trials, n) - 0.5;
uTie = rand(trials, 1);
echo = exp(1i * phase) .* s;
hits = false(numel(snrDb), trials);
for j = 1:numel(snrDb)
    amp = 10^(snrDb(j) / 20);
    z = radarwave.applyReceiver(rx, re + amp * real(echo), im + amp * imag(echo), uRe, uIm);
    stat = abs(z * s').^2;           % s' is the conjugate transpose: sum z .* conj(s)
    hits(j, :) = radarwave.decide(thr, stat, uTie).';
end
end
