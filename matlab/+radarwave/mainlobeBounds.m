function [left, right] = mainlobeBounds(mag, p)
%MAINLOBEBOUNDS Indices of the first minimum on each side of the peak p (inclusive).
mag = abs(mag(:));
if nargin < 2
    [~, p] = max(mag);
end
tol = 1e-12 * mag(p);
left = p;
while left > 1 && mag(left - 1) < mag(left) + tol
    left = left - 1;
end
right = p;
while right < numel(mag) && mag(right + 1) < mag(right) + tol
    right = right + 1;
end
end
