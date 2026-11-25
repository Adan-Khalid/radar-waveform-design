function z = applyReceiver(rx, re, im, uRe, uIm)
%APPLYRECEIVER Quantize I and Q (real arrays) with the receiver's ADC; complex output.
%   Noise has power 1 (each component standard deviation 1/sqrt(2)).
if isinf(rx.bits)
    z = re + 1i * im;
    return
end
A = rx.aOverSigma / sqrt(2);
z = radarwave.quantizeDither(re, rx.bits, A, rx.kind, rx.dither, uRe) + ...
    1i * radarwave.quantizeDither(im, rx.bits, A, rx.kind, rx.dither, uIm);
end
