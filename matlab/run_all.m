function run_all(skipAdc)
%RUN_ALL Run the MATLAB tests, then every study, one after another.
%   How to run it in MATLAB Online (free licence is enough):
%   1. Upload to MATLAB Drive, keeping the folder names: matlab/ and reports/
%      (reports/ holds the Python results and reports/crosscheck/ the exported values).
%   2. In MATLAB Online:  cd radar-waveform-design/matlab   then   run_all
%   3. Results and figures are written to reports/matlab/. Download that folder.
%   run_all(true) skips the Monte Carlo ADC study (the slowest part, a few minutes).
%   One study that fails prints its error and the others still run.
if nargin < 1
    skipAdc = false;
end
here = fileparts(mfilename('fullpath'));
addpath(here);                       % makes the +radarwave package visible
addpath(fullfile(here, 'studies'));

fprintf('\n=== Tests ===\n');
try
    results = runtests(fullfile(here, 'tests'));
    disp(table(results));
    fprintf('%d passed, %d failed, %d incomplete\n', sum([results.Passed]), ...
        sum([results.Failed]), sum([results.Incomplete]));
catch err
    fprintf('Tests could not run: %s\n', err.message);
end

studies = {'study_waveforms', 'study_ambiguity', 'study_sidelobes', 'study_adc', ...
    'study_crosscheck'};
if skipAdc
    studies(strcmp(studies, 'study_adc')) = [];
end
for i = 1:numel(studies)
    fprintf('\n=== %s ===\n', studies{i});
    t0 = tic;
    try
        runOne(fullfile(here, 'studies', [studies{i} '.m']));
        fprintf('(%s took %.0f s)\n', studies{i}, toc(t0));
    catch err
        fprintf('FAILED: %s\n', err.message);
        if ~isempty(err.stack)
            fprintf('  at %s, line %d\n', err.stack(1).name, err.stack(1).line);
        end
    end
end
end

function runOne(scriptFile)
% Each study runs in its own workspace, so the studies cannot disturb each other.
run(scriptFile);
end
