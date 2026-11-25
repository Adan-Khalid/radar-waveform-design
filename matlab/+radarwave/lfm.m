function s = lfm(timeBandwidth, samplesPerCell)
%LFM Linear FM pulse sweeping -B/2 to +B/2, unit energy (column vector).
%   s = radarwave.lfm(BT, samplesPerCell)
%   s(t) = exp(j*pi*(B/T)*(t - T/2)^2), 0 <= t < T, sampled at the centre of each of
%   round(BT*samplesPerCell) samples. Same as radarwave.waveforms.lfm in Python.
if nargin < 2
    samplesPerCell = 1;
end
n = round(timeBandwidth * samplesPerCell);
t = ((0:n-1).' + 0.5) / n;   % fraction of the pulse length
s = radarwave.unitEnergy(exp(1i * pi * timeBandwidth * (t - 0.5).^2));
end
