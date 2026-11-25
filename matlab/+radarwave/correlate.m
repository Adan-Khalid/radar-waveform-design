function [y, lag] = correlate(x, ref, mode)
%CORRELATE Cross-correlation by FFT: y(k) = sum_n x(n + k) * conj(ref(n)).
%   [y, lag] = radarwave.correlate(x, ref)          lags -(numel(ref)-1) .. numel(x)-1
%   [y, lag] = radarwave.correlate(X, ref, 'valid')  lags 0 .. size(X,2)-numel(ref)
%   An echo of ref starting at sample d of x gives the matched peak at lag d.
%   x may be a vector, or a matrix with one record per ROW (same as Python's last axis).
if nargin < 3
    mode = 'full';
end
ref = ref(:).';
isVec = isvector(x);
if isVec
    x = x(:).';
end
m = size(x, 2);
nr = numel(ref);
nfft = 2^nextpow2(m + nr - 1);
c = ifft(fft(x, nfft, 2) .* conj(fft(ref, nfft, 2)), [], 2);
switch mode
    case 'full'
        y = [c(:, nfft-nr+2:nfft), c(:, 1:m)];
        lag = -(nr-1):(m-1);
    case 'valid'
        if m < nr
            error('radarwave:correlate', 'x is shorter than the reference.');
        end
        y = c(:, 1:m-nr+1);
        lag = 0:(m-nr);
    otherwise
        error('radarwave:correlate', 'mode must be full or valid.');
end
if isVec
    y = y(:);
    lag = lag(:);
end
end
