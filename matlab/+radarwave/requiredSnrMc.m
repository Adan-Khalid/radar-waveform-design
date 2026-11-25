function curve = requiredSnrMc(s, rx, thr, opts, predictedDb)
%REQUIREDSNRMC Coarse grid with few trials to bracket Pd = target, then a fine grid
%   with all trials. Returns struct snrDb (column), hits (numel(snrDb) x trials).
%   opts fields: pdTarget, trials, coarseTrials, seed, fineHalfWidth, fineStep,
%   maxSnrDb (see study_adc.m).
lo = predictedDb - 3;
hi = min(predictedDb + 9, opts.maxSnrDb);
approx = NaN;
while isnan(approx)
    grid = (lo:0.5:hi).';
    pd = mean(radarwave.detectionsAt(s, rx, thr, grid, opts.coarseTrials, opts.seed + 1), 2);
    x = radarwave.crossingDb(grid, pd, opts.pdTarget);
    if isfinite(x)
        approx = x;
    elseif x == -Inf
        hi = lo + 0.5;
        lo = lo - 6;
    elseif hi >= opts.maxSnrDb
        approx = Inf;
    else
        lo = hi - 0.5;
        hi = min(hi + 9, opts.maxSnrDb);
    end
end
if ~isfinite(approx)
    grid = (opts.maxSnrDb - 2 * opts.fineHalfWidth:opts.fineStep:opts.maxSnrDb).';
else
    grid = round((approx - opts.fineHalfWidth:opts.fineStep:approx + opts.fineHalfWidth).' * 1e6) / 1e6;
end
hits = radarwave.detectionsAt(s, rx, thr, grid, opts.trials, opts.seed);
curve = struct('snrDb', grid, 'hits', hits);
end
