function v = paprDb(s)
%PAPRDB Peak-to-average power ratio of the samples, dB.
p = abs(s(:)).^2;
v = 10 * log10(max(p) / mean(p));
end
