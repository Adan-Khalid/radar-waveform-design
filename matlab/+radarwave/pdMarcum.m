function pd = pdMarcum(snrDb, pfa)
%PDMARCUM Pd of a non-fluctuating target with random phase, square-law detector, one
%   pulse: Pd = Q1(sqrt(2*SNR), sqrt(-2*ln(Pfa))) (Marcum's Q; marcumq is in the
%   Signal Processing Toolbox).
snr = 10.^(snrDb / 10);
pd = marcumq(sqrt(2 * snr), sqrt(-2 * log(pfa)));
end
