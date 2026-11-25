function sym = ofdmSymbols(n)
%OFDMSYMBOLS The QPSK symbols Python used for the n-subcarrier OFDM pulse.
file = fullfile(radarwave.repoRoot(), 'reports', 'crosscheck', 'ofdm_symbols.csv');
if ~isfile(file)
    error('radarwave:ofdm', ['Missing %s. Run "python scripts/export_crosscheck.py" ' ...
        'and upload reports/crosscheck with the matlab folder.'], file);
end
T = readtable(file);
sel = T.subcarriers == n;
if ~any(sel)
    error('radarwave:ofdm', 'No exported symbols for %d subcarriers.', n);
end
sym = complex(T.re(sel), T.im(sel));
end
