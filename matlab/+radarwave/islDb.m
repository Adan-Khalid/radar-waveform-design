function v = islDb(y)
%ISLDB Integrated sidelobe level: sidelobe energy over mainlobe energy, dB.
mag = abs(y(:));
[~, p] = max(mag);
[l, r] = radarwave.mainlobeBounds(mag, p);
side = [mag(1:l-1); mag(r+1:end)];
v = 10 * log10(sum(side.^2) / sum(mag(l:r).^2));
end
