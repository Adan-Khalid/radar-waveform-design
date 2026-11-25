function [x, w] = gaussLegendre(n)
%GAUSSLEGENDRE Nodes and weights on [-1, 1] (Golub-Welsch), base MATLAB only.
k = 1:n-1;
beta = k ./ sqrt(4 * k.^2 - 1);
[V, D] = eig(diag(beta, 1) + diag(beta, -1));
[x, order] = sort(diag(D));
w = 2 * V(1, order).'.^2;
end
