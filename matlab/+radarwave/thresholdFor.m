function thr = thresholdFor(stats, pfa)
%THRESHOLDFOR Neyman-Pearson threshold, randomised on ties (struct value, gamma).
%   Decide "target" when stat > value, and with probability gamma when stat equals
%   value (within a relative 1e-4), so the false-alarm rate is pfa even when the
%   statistic takes only a few values (one bit, short binary code).
tol = 1e-4;
v = sort(stats(:));
n = numel(v);
t = v(ceil((n - 1) * (1 - pfa)) + 1);   % an observed value ("higher" quantile)
band = tol * t;
nAbove = sum(v > t + band);
nTie = sum(abs(v - t) <= band);
if nTie > 0
    gamma = min(max((pfa * n - nAbove) / nTie, 0), 1);
else
    gamma = 0;
end
thr = struct('value', t, 'gamma', gamma);
end
