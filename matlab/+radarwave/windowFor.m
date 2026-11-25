function w = windowFor(name, n, sidelobeDb, nbar)
%WINDOWFOR Symmetric amplitude window by name (column vector):
%   'rectangular', 'hann', 'hamming', 'blackman', 'taylor' (sidelobeDb, nbar).
if nargin < 3, sidelobeDb = 35; end
if nargin < 4, nbar = 5; end
i = (0:n-1).';
c = @(k) cos(2 * pi * k * i / max(n - 1, 1));
switch lower(name)
    case 'rectangular'
        w = ones(n, 1);
    case 'hann'
        w = 0.5 - 0.5 * c(1);
    case 'hamming'
        w = 0.54 - 0.46 * c(1);
    case 'blackman'
        w = 0.42 - 0.5 * c(1) + 0.08 * c(2);
    case 'taylor'
        w = radarwave.taylorWindow(n, sidelobeDb, nbar);
    otherwise
        error('radarwave:window', 'Unknown window %s.', name);
end
end
