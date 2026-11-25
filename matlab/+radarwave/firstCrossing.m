function v = firstCrossing(x, y, level)
%FIRSTCROSSING Smallest x where y first falls below level (linear interpolation).
j = find(y < level, 1);
if isempty(j)
    v = Inf;
elseif j == 1
    v = x(1);
else
    v = x(j-1) + (level - y(j-1)) * (x(j) - x(j-1)) / (y(j) - y(j-1));
end
end
