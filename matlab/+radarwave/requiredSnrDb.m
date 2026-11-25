function snrDb = requiredSnrDb(pd, pfa)
%REQUIREDSNRDB SNR (dB) at which radarwave.pdMarcum reaches pd.
opts = optimset('TolX', 1e-10);
snrDb = fzero(@(x) radarwave.pdMarcum(x, pfa) - pd, [-20, 40], opts);
end
