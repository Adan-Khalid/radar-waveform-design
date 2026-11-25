function s = unitEnergy(s)
%UNITENERGY Scale a waveform to unit energy, sum(|s|^2) = 1 (column vector).
s = s(:);
s = s / sqrt(sum(abs(s).^2));
end
