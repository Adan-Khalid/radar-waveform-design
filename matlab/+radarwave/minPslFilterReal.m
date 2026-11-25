function ref = minPslFilterReal(s, len)
%MINPSLFILTERREAL Minimum peak sidelobe filter for a real code (linear programme).
%   minimise t subject to |z_k| <= t (k ~= c), z_c = 1, z = S*g, g real
%   (Baden and Cohen). Needs linprog (Optimization Toolbox, in free MATLAB Online).
s = real(s(:));
S = radarwave.convolutionMatrix(s, len);
nOut = size(S, 1);
c = floor((nOut - 1) / 2) + 1;
rest = S([1:c-1, c+1:nOut], :);
m = size(rest, 1);
f = [zeros(len, 1); 1];
A = [rest, -ones(m, 1); -rest, -ones(m, 1)];
b = zeros(2 * m, 1);
Aeq = [S(c, :), 0];
beq = 1;
lb = [-Inf(len, 1); 0];
ub = Inf(len + 1, 1);
opts = optimoptions('linprog', 'Display', 'none');
x = linprog(f, A, b, Aeq, beq, lb, ub, opts);
ref = flipud(x(1:len));
end
