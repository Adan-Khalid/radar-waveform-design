function label = waveformLabel(name, timeBandwidth)
%WAVEFORMLABEL Display name used in every table, e.g. "P4 64" or "LFM (BT 64)".
switch lower(name)
    case 'lfm',    label = sprintf("LFM (BT %d)", timeBandwidth);
    case 'barker', label = sprintf("Barker %d", timeBandwidth);
    case 'frank',  label = sprintf("Frank %d", timeBandwidth);
    case 'p4',     label = sprintf("P4 %d", timeBandwidth);
    case 'ofdm',   label = sprintf("OFDM %d", timeBandwidth);
end
end
