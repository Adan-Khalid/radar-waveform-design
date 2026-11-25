function ph = frankPhases(n)
%FRANKPHASES Frank code of length n = M^2: phase(i,k) = 2*pi*i*k/M, chip i*M + k.
m = round(sqrt(n));
if m * m ~= n
    error('radarwave:frank', 'A Frank code needs a square length.');
end
[K, I] = meshgrid(0:m-1, 0:m-1);   % I(r,c) = r-1, K(r,c) = c-1
P = 2 * pi * I .* K / m;
ph = reshape(P.', [], 1);          % row by row: chip index i*M + k
end
