function s = ofdm(symbols, samplesPerCell)
%OFDM OFDM pulse: N subcarriers spaced 1/T, centred on zero, no cyclic prefix.
%   s(t) = sum_k X_k exp(j*2*pi*(k - (N-1)/2)*t/T), so B = N/T and BT = N.
%   Use the symbols exported by Python (reports/crosscheck/ofdm_symbols.csv) to get
%   exactly the same pulse, since MATLAB and NumPy random numbers differ.
if nargin < 2
    samplesPerCell = 1;
end
symbols = symbols(:);
n = numel(symbols);
t = ((0:n*samplesPerCell-1).' + 0.5) / (n * samplesPerCell);
f = (0:n-1) - (n - 1) / 2;
s = radarwave.unitEnergy(exp(2i * pi * t * f) * symbols);
end
