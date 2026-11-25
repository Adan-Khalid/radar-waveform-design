function rx = makeReceiver(bits, kind, dither, aOverSigma)
%MAKERECEIVER ADC settings of one receiver (struct).
%   rx = radarwave.makeReceiver()                 unquantized (bits = Inf)
%   rx = radarwave.makeReceiver(3)                3 bits, mid-rise, no dither, Max's step
%   rx = radarwave.makeReceiver(3, 'midrise', 'subtractive', aOverSigma)
%   aOverSigma: full scale over the noise standard deviation of one component
%   (default: Max's Gaussian optimum for the bit depth).
if nargin < 1, bits = Inf; end
if nargin < 2 || isempty(kind), kind = 'midrise'; end
if nargin < 3 || isempty(dither), dither = 'none'; end
if nargin < 4 || isempty(aOverSigma)
    if isinf(bits)
        aOverSigma = NaN;
    else
        aOverSigma = radarwave.gaussianOptimalFullScale(bits);
    end
end
if isinf(bits)
    label = "unquantized";
else
    label = sprintf("%d bit", bits);
    if ~strcmp(kind, 'midrise'), label = label + ", " + kind; end
    if ~strcmp(dither, 'none'), label = label + ", " + dither + " dither"; end
end
rx = struct('bits', bits, 'kind', kind, 'dither', dither, 'aOverSigma', aOverSigma, ...
    'label', label);
end
