function s = phaseCode(phases, samplesPerChip)
%PHASECODE Rectangular chips with the given phases, unit energy (column vector).
if nargin < 2
    samplesPerChip = 1;
end
chips = exp(1i * phases(:));
s = radarwave.unitEnergy(repelem(chips, samplesPerChip));
end
