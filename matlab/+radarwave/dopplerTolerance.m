function [summ, cut] = dopplerTolerance(s, ref, nuTMax, step, samplesPerCell)
%DOPPLERTOLERANCE Doppler tolerance numbers, as in Python's doppler_tolerance:
%   nu_t_1db, nu_t_3db (first crossing), loss at nu*T = 0.5, 1, 2, 4 (dB),
%   peak shift at nu*T = 1 and 4 (cells), PSL at nu*T = 0.25 and 0.5.
if nargin < 2, ref = []; end
if nargin < 3, nuTMax = 4; end
if nargin < 4, step = 0.01; end
if nargin < 5, samplesPerCell = 1; end
nu = round((0:step:nuTMax + step/2) * 1e10) / 1e10;
nu = nu(nu <= nuTMax + step/2);
cut = radarwave.dopplerCutSummary(s, ref, nu, samplesPerCell);
at = @(field, v) pick(nu, cut.(field), v, step);
summ = struct();
summ.nu_t_1db = radarwave.firstCrossing(nu, cut.peak_loss_db, -1);
summ.nu_t_3db = radarwave.firstCrossing(nu, cut.peak_loss_db, -3);
summ.loss_at_0p5_db = at('peak_loss_db', 0.5);
summ.loss_at_1_db = at('peak_loss_db', 1);
summ.loss_at_2_db = at('peak_loss_db', 2);
summ.loss_at_4_db = at('peak_loss_db', 4);
summ.shift_at_1_cells = at('shift_cells', 1);
summ.shift_at_4_cells = at('shift_cells', 4);
summ.psl_at_0p25_db = at('psl_db', 0.25);
summ.psl_at_0p5_db = at('psl_db', 0.5);
end

function v = pick(nu, values, target, step)
[d, j] = min(abs(nu - target));
if d < step / 2
    v = values(j);
else
    v = NaN;
end
end
