function v = pslDb(y)
%PSLDB Peak sidelobe level: largest value outside the mainlobe, dB below the peak.
mag = abs(y(:));
[pk, p] = max(mag);
[l, r] = radarwave.mainlobeBounds(mag, p);
side = [mag(1:l-1); mag(r+1:end)];
if isempty(side) || max(side) == 0
    v = -Inf;
else
    v = 20 * log10(max(side) / pk);
end
end
