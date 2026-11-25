function [lossTable, pfaTable] = runReceivers(s, receivers, opts)
%RUNRECEIVERS Thresholds, Pfa check, required SNR and paired losses for receivers.
%   The first receiver is the reference for the losses (normally unquantized).
%   opts fields: pfa, pdTarget, noiseStats, trials, coarseTrials, bootstrap, seed,
%   fineHalfWidth, fineStep, maxSnrDb.
nRx = numel(receivers);
stats = radarwave.noiseStatistics(s, receivers, opts.noiseStats, opts.seed + 11);
fresh = radarwave.noiseStatistics(s, receivers, round(opts.noiseStats / 2), opts.seed + 12);
rng(opts.seed + 13, 'twister');
thr = cell(nRx, 1);
pfaRows = cell(nRx, 7);
for k = 1:nRx
    thr{k} = radarwave.thresholdFor(stats{k}, opts.pfa);
    u = rand(numel(fresh{k}), 1);
    fa = sum(radarwave.decide(thr{k}, fresh{k}, u));
    [plo, phi] = radarwave.wilson(fa, numel(fresh{k}));
    pfaRows(k, :) = {receivers(k).label, fa, numel(fresh{k}), fa / numel(fresh{k}), plo, phi, ...
        plo <= opts.pfa && opts.pfa <= phi};
end
pfaTable = cell2table(pfaRows, 'VariableNames', ...
    {'receiver', 'false_alarms', 'statistics', 'pfa_measured', 'pfa_lo', 'pfa_hi', 'design_inside'});

req0 = radarwave.requiredSnrDb(opts.pdTarget, opts.pfa);
curves = cell(nRx, 1);
for k = 1:nRx
    rx = receivers(k);
    if isinf(rx.bits)
        pred = req0;
    else
        eff = radarwave.quantizerEfficiency(rx.bits, rx.aOverSigma, rx.kind, rx.dither);
        pred = req0 + min(-10 * log10(max(eff, 1e-12)), 40);
    end
    curves{k} = radarwave.requiredSnrMc(s, rx, thr{k}, opts, pred);
    fprintf('  %-34s done\n', rx.label);
end

% Paired bootstrap: Poisson(1) weights on trials, the same for every receiver.
rng(opts.seed + 14, 'twister');
W = poissrnd(1, opts.trials, opts.bootstrap);
boot = zeros(nRx, opts.bootstrap);
req = zeros(nRx, 1);
for k = 1:nRx
    c = curves{k};
    req(k) = radarwave.crossingDb(c.snrDb, mean(c.hits, 2), opts.pdTarget);
    pdRep = (double(c.hits) * W) ./ sum(W, 1);
    boot(k, :) = radarwave.crossingDb(c.snrDb, pdRep, opts.pdTarget);
end
rows = cell(nRx, 9);
for k = 1:nRx
    rx = receivers(k);
    diffs = boot(k, :) - boot(1, :);
    ci = prctile(boot(k, :), [2.5, 97.5]);
    lci = prctile(diffs, [2.5, 97.5]);
    if isinf(rx.bits)
        theoryLoss = 0;
    else
        theoryLoss = -10 * log10(radarwave.quantizerEfficiency(rx.bits, rx.aOverSigma, rx.kind, rx.dither));
    end
    rows(k, :) = {rx.label, req(k), ci(1), ci(2), req(k) - req(1), lci(1), lci(2), theoryLoss, ...
        thr{k}.value};
end
lossTable = cell2table(rows, 'VariableNames', {'receiver', 'required_snr_db', ...
    'required_snr_lo', 'required_snr_hi', 'loss_db', 'loss_lo', 'loss_hi', ...
    'theory_loss_db', 'threshold'});
end
