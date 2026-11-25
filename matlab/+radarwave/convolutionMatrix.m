function S = convolutionMatrix(s, len)
%CONVOLUTIONMATRIX S such that S*g is the full convolution of s with g (len taps).
s = s(:);
n = numel(s);
nOut = n + len - 1;
I = (0:nOut-1).' - (0:len-1);
S = zeros(nOut, len);
v = I >= 0 & I < n;
S(v) = s(I(v) + 1);
end
