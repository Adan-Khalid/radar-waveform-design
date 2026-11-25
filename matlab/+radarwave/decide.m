function d = decide(thr, stat, u)
%DECIDE Detections from statistics with a randomised threshold (see thresholdFor).
band = 1e-4 * thr.value;
d = stat > thr.value + band;
if thr.gamma > 0
    d = d | (abs(stat - thr.value) <= band & u < thr.gamma);
end
end
