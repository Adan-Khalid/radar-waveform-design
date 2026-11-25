function s = makeWaveform(name, timeBandwidth, samplesPerCell)
%MAKEWAVEFORM Waveform of one design point by name: 'lfm', 'barker', 'frank', 'p4',
%   'ofdm'. Code length and subcarrier count equal the time-bandwidth product.
%   OFDM symbols are read from reports/crosscheck/ofdm_symbols.csv (written by
%   scripts/export_crosscheck.py), so MATLAB builds exactly the pulses Python uses.
if nargin < 3, samplesPerCell = 1; end
switch lower(name)
    case 'lfm'
        s = radarwave.lfm(timeBandwidth, samplesPerCell);
    case 'barker'
        s = radarwave.phaseCode(radarwave.barkerPhases(timeBandwidth), samplesPerCell);
    case 'frank'
        s = radarwave.phaseCode(radarwave.frankPhases(timeBandwidth), samplesPerCell);
    case 'p4'
        s = radarwave.phaseCode(radarwave.p4Phases(timeBandwidth), samplesPerCell);
    case 'ofdm'
        s = radarwave.ofdm(radarwave.ofdmSymbols(timeBandwidth), samplesPerCell);
    otherwise
        error('radarwave:waveform', 'Unknown waveform %s.', name);
end
end
