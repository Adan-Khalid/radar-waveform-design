function x = crossingDb(snrDb, pd, target)
%CROSSINGDB SNR where Pd reaches target: linear interpolation above the last point
%   below it. pd may have one column per bootstrap replicate. Inf if never reached,
%   -Inf if the first point is already above the target.
if isvector(pd), pd = pd(:); end
snrDb = snrDb(:);
nRep = size(pd, 2);
x = zeros(1, nRep);
for r = 1:nRep
    j = find(pd(:, r) < target, 1, 'last');
    if isempty(j)
        x(r) = -Inf;
    elseif j == numel(snrDb)
        x(r) = Inf;
    else
        p0 = pd(j, r);
        p1 = pd(j + 1, r);
        x(r) = snrDb(j) + (target - p0) * (snrDb(j + 1) - snrDb(j)) / (p1 - p0);
    end
end
end
