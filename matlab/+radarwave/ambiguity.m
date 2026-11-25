function [chi, lags, nuT] = ambiguity(s, ref, nuT, lags)
%AMBIGUITY Narrowband (cross-)ambiguity function from the definition.
%   chi(m, k) = sum_n s(n) exp(j*2*pi*nuT(k)*n/N) conj(ref(n - m)),  N = numel(s)
%   Rows are lags m (default -(numel(ref)-1) .. numel(s)-1), columns Doppler values
%   given as nu*T (cycles over the pulse). ref = [] gives the auto-ambiguity function.
%   With nuT = (0:K-1)*N/K, K >= N (a full DFT grid), sum(|chi|^2)/K = Es*Eref exactly.
s = s(:);
if nargin < 2 || isempty(ref), ref = s; end
ref = ref(:);
ns = numel(s);
nr = numel(ref);
if nargin < 3 || isempty(nuT), nuT = 0:ns-1; end
if nargin < 4 || isempty(lags), lags = -(nr-1):(ns-1); end
lags = lags(:);
n = 0:ns-1;
idx = n - lags;                         % implicit expansion: (numLags x ns)
valid = idx >= 0 & idx < nr;
R = zeros(size(idx));
R(valid) = conj(ref(idx(valid) + 1));
P = R .* s.';
D = exp(2i * pi * (n.') * nuT(:).' / ns);
chi = P * D;
end
