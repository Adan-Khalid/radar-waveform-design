function w = taylorWindow(n, sidelobeDb, nbar)
%TAYLORWINDOW Taylor window written out from its formula (same as taylorwin and
%   SciPy's taylor with norm=False). sidelobeDb is positive, e.g. 35.
if nargin < 2, sidelobeDb = 35; end
if nargin < 3, nbar = 5; end
A = acosh(10^(sidelobeDb / 20)) / pi;
sp2 = nbar^2 / (A^2 + (nbar - 0.5)^2);
m = (1:nbar-1).';
F = zeros(nbar - 1, 1);
for j = 1:nbar-1
    mj = m(j);
    numer = prod(1 - mj^2 / sp2 ./ (A^2 + (m - 0.5).^2));
    others = m(m ~= mj);
    denom = 2 * prod(1 - mj^2 ./ others.^2);
    F(j) = (-1)^(mj + 1) * numer / denom;
end
x = ((0:n-1).' - n/2 + 0.5) / n;
w = 1 + 2 * cos(2 * pi * x * m.') * F;
end
