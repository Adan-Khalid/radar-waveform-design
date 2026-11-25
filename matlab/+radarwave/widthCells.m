function w = widthCells(y, levelDb, samplesPerCell)
%WIDTHCELLS Mainlobe width at levelDb below the peak, in range cells (units of 1/B).
%   Crossings by linear interpolation of |y| between samples.
if nargin < 2, levelDb = -3; end
if nargin < 3, samplesPerCell = 1; end
mag = abs(y(:));
[pk, p] = max(mag);
level = pk * 10^(levelDb / 20);
w = (crossing(mag, p, level, 1) - crossing(mag, p, level, -1)) / samplesPerCell;
end

function x = crossing(mag, p, level, direction)
i = p;
while i + direction >= 1 && i + direction <= numel(mag) && mag(i + direction) > level
    i = i + direction;
end
j = i + direction;
if j < 1 || j > numel(mag)
    x = i;
    return
end
frac = (mag(i) - level) / (mag(i) - mag(j));
x = i + direction * frac;
end
