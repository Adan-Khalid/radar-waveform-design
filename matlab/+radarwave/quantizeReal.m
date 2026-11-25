function q = quantizeReal(x, bits, fullScale, kind)
%QUANTIZEREAL Uniform quantizer on real samples. bits = Inf returns x unchanged.
%   Step d = 2*A/2^bits for both kinds:
%   'midrise'  (default): 2^bits levels at +-d/2, +-3d/2, ... (one bit = sign)
%   'midtread': 2^bits - 1 levels at 0, +-d, ... (needs 2 or more bits)
if nargin < 4, kind = 'midrise'; end
if isinf(bits)
    q = x;
    return
end
d = 2 * fullScale / 2^bits;
switch kind
    case 'midrise'
        half = 2^(bits - 1);
        k = min(max(floor(x / d), -half), half - 1);
        q = (k + 0.5) * d;
    case 'midtread'
        if bits < 2
            error('radarwave:quantize', 'A mid-tread quantizer needs at least 2 bits.');
        end
        kmax = 2^(bits - 1) - 1;
        k = min(max(floor(x / d + 0.5), -kmax), kmax);
        q = k * d;
    otherwise
        error('radarwave:quantize', 'kind must be midrise or midtread.');
end
end
