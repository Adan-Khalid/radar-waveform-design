function A = gaussianOptimalFullScale(bits, sigma)
%GAUSSIANOPTIMALFULLSCALE Full scale giving Max's (1960, Table II) optimal mid-rise
%   step for a Gaussian input of standard deviation sigma.
if nargin < 2, sigma = 1; end
steps = [1.596, 0.9957, 0.5860, 0.3352];   % 1, 2, 3, 4 bits, unit variance
A = 2^(bits - 1) * steps(bits) * sigma;
end
