function v = mismatchLossDb(s, ref)
%MISMATCHLOSSDB SNR loss of a receive reference against the matched filter, dB.
%   v = -10*log10( max|y|^2 / (||s||^2 ||ref||^2) ), y = correlate(s, ref). Never negative.
y = radarwave.correlate(s, ref);
v = -10 * log10(max(abs(y))^2 / (sum(abs(s).^2) * sum(abs(ref).^2)));
end
