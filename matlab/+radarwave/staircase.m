function [thresholds, levels] = staircase(bits, fullScale, kind)
%STAIRCASE Decision thresholds (ascending, row) and output levels (row) of the quantizer.
if nargin < 3, kind = 'midrise'; end
d = 2 * fullScale / 2^bits;
switch kind
    case 'midrise'
        half = 2^(bits - 1);
        levels = ((-half:half-1) + 0.5) * d;
        thresholds = (-half+1:half-1) * d;
    case 'midtread'
        kmax = 2^(bits - 1) - 1;
        levels = (-kmax:kmax) * d;
        thresholds = ((-kmax:kmax-1) + 0.5) * d;
    otherwise
        error('radarwave:quantize', 'kind must be midrise or midtread.');
end
end
