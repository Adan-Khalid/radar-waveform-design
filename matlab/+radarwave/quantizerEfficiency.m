function eta = quantizerEfficiency(bits, aOverSigma, kind, dither, nodes)
%QUANTIZEREFFICIENCY Low-SNR efficiency of a quantizer on a Gaussian input (Bussgang):
%   eta = E[x z]^2 / (sigma^2 E[z^2]), z the receiver output after dither (and after
%   removing it for subtractive dither). Predicted loss: -10*log10(eta).
%   aOverSigma: full scale A over the noise standard deviation of one component.
%   Expectations over x are exact (normal CDF sums); the dither average uses
%   Gauss-Legendre nodes. Same as Python's theory.quantizer_efficiency.
if nargin < 3, kind = 'midrise'; end
if nargin < 4, dither = 'none'; end
if nargin < 5, nodes = 64; end
if isinf(bits)
    eta = 1;
    return
end
[thr, lev] = radarwave.staircase(bits, aOverSigma, kind);
if strcmp(dither, 'none')
    [exq, ~, eq2] = moments(thr, lev, 0);
    if eq2 <= 0
        eta = 0;
    else
        eta = exq^2 / eq2;
    end
    return
end
d = 2 * aOverSigma / 2^bits;
[g, w] = radarwave.gaussLegendre(nodes);
shifts = g * d / 2;
w = w / 2;
[exq, eq, eq2] = moments(thr, lev, shifts);
exz = sum(w .* exq);
switch dither
    case 'nonsubtractive'
        ez2 = sum(w .* eq2);
    case 'subtractive'
        ez2 = sum(w .* (eq2 - 2 * shifts .* eq + shifts.^2));
    otherwise
        error('radarwave:dither', 'dither must be none, nonsubtractive or subtractive.');
end
eta = exz^2 / ez2;
end

function [exq, eq, eq2] = moments(thr, lev, shift)
% E[x q(x+shift)], E[q(x+shift)], E[q(x+shift)^2] for x ~ N(0,1); shift a column.
shift = shift(:);
lo = [-Inf, thr] - shift;
hi = [thr, Inf] - shift;
Phi = @(z) 0.5 * erfc(-z / sqrt(2));
phi = @(z) exp(-0.5 * z.^2) / sqrt(2 * pi);
prob = Phi(hi) - Phi(lo);
exq = sum(lev .* (phi(lo) - phi(hi)), 2);
eq = sum(lev .* prob, 2);
eq2 = sum(lev.^2 .* prob, 2);
end
