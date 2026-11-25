function ref = lsMismatchedFilter(s, len, loading)
%LSMISMATCHEDFILTER Least-squares mismatched filter, returned as a receive reference.
%   Output z = S*g. Fix z(c) = 1 at the centre and minimise the energy of the rest:
%   g = R \ a' / (a * (R \ a')),  R = Srest'*Srest + loading*I  (Ackroyd and Ghani).
%   radarwave.correlate(x, ref) then equals the convolution of x with g.
if nargin < 3, loading = 1e-9; end
s = s(:);
S = radarwave.convolutionMatrix(s, len);
nOut = size(S, 1);
c = floor((nOut - 1) / 2) + 1;
a = S(c, :);
rest = S([1:c-1, c+1:nOut], :);
R = rest' * rest + loading * sum(abs(s).^2) * eye(len);
g = R \ a';
g = g / (a * g);
ref = conj(flipud(g));
end
