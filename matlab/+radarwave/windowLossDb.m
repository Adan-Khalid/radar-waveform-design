function v = windowLossDb(w)
%WINDOWLOSSDB SNR loss of an amplitude weighting on a constant-envelope pulse, dB:
%   -10*log10( mean(w)^2 / mean(w.^2) ).
w = w(:);
v = -10 * log10(mean(w)^2 / mean(w.^2));
end
