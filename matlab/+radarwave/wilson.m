function [lo, hi] = wilson(k, n, z)
%WILSON Wilson score interval for a binomial proportion k/n (95 % by default).
if nargin < 3, z = 1.959964; end
p = k ./ n;
den = 1 + z^2 ./ n;
centre = (p + z^2 ./ (2 * n)) ./ den;
half = z * sqrt(p .* (1 - p) ./ n + z^2 ./ (4 * n.^2)) ./ den;
lo = centre - half;
hi = centre + half;
end
