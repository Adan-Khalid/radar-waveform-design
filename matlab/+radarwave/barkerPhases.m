function ph = barkerPhases(n)
%BARKERPHASES Barker code of length n as phases 0 or pi (column vector).
switch n
    case 2,  c = [1 -1];
    case 3,  c = [1 1 -1];
    case 4,  c = [1 1 -1 1];
    case 5,  c = [1 1 1 -1 1];
    case 7,  c = [1 1 1 -1 -1 1 -1];
    case 11, c = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
    case 13, c = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    otherwise
        error('radarwave:barker', 'No Barker code of length %d.', n);
end
ph = pi * (c(:) < 0);
end
