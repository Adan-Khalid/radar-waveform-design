function q = quantizeDither(x, bits, fullScale, kind, dither, u)
%QUANTIZEDITHER Quantize real samples with uniform dither of one step peak to peak.
%   u: uniform values in [-1/2, 1/2), same size as x.
%   dither 'none': q(x); 'nonsubtractive': q(x + d); 'subtractive': q(x + d) - d.
if isinf(bits) || strcmp(dither, 'none')
    q = radarwave.quantizeReal(x, bits, fullScale, kind);
    return
end
dv = u * (2 * fullScale / 2^bits);
q = radarwave.quantizeReal(x + dv, bits, fullScale, kind);
switch dither
    case 'nonsubtractive'
        % keep q
    case 'subtractive'
        q = q - dv;
    otherwise
        error('radarwave:dither', 'dither must be none, nonsubtractive or subtractive.');
end
end
