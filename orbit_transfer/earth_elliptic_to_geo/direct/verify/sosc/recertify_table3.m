function T = recertify_table3(thrustList)
% RECERTIFY_TABLE3  Batch re-certify the campaign's headline Table-3 rows with
% the SOSC certificate (FINAL direct reduced-Hessian method, DESIGN sec 12).
% Writes a SIDECAR verdict file per rung (results/sosc/sosc_<tag>.mat, variable
% `sosc`) and prints a summary table; the campaign `results/*.mat` caches are
% never written to (read-only here).
%
% INPUTS:
%   thrustList - vector of thrust levels [N] (e.g. [10 5 2.5 1 0.5])
%
% OUTPUTS:
%   T - struct array, one row per certified thrust level:
%       .thrustN    - thrust level [N]
%       .tag        - campaign row tag [char]
%       .verdict    - one of PASS/WEAK_MIN/FAIL/INCONCLUSIVE/ERROR [char]
%       .status     - tiered-gate status string [char]
%       .nFlat      - dimension of the flat manifold (reduced nzero) [int]
%       .sensStable - reduced negative-count zt-stable across the band [bool]
%       .method     - inertia method, 'reduced-eig' or 'scale-skip' [char]
%       .robust     - true iff the reduced inertia was computed [bool]
%       .red        - reduced-inertia struct .npos .nneg .nzero
%       .redMinEig  - min eigenvalue of the reduced Hessian [double]
%       .drift      - warm re-solve drift [double]
%       .stat       - KKT stationarity residual (sosc.kkt.stat) [double]
%     Rows for MISSING campaign files are skipped (not appended to T).
%
% REFERENCES: process/DESIGN_sosc.md sec 8 (plug-in points), sec 12
%   (FINAL inertia method + verdict taxonomy), sec 11.6 (tiered gate).
resDir  = fullfile(module_root(),'results');
sideDir = fullfile(resDir,'sosc');
if ~isfolder(sideDir), mkdir(sideDir); end

% tagOf: the CERTIFIED HEADLINE row per rung. 10/5/2.5 N are the MEE_M2 fuel
% rows; 1 N and 0.5 N headline numbers (1371.44 kg, 1375.28 kg) are the
% PSR-refined solutions, so certify those PSR-final rows (sosc_load_row
% normalizes both the `res` and `out` shapes).
tagOf = containers.Map({10,5,2.5,1,0.5}, ...
    {'MEE_M2_10N','MEE_M2_5N','MEE_M2_2p5N', ...
     'MEE_M2_1N_PSR_psr_final','MEE_M2_0p5N_PSR_psr_final'});

T = struct('thrustN',{},'tag',{},'verdict',{},'status',{},'nFlat',{}, ...
    'sensStable',{},'method',{},'robust',{},'red',{},'redMinEig',{}, ...
    'drift',{},'stat',{});

fprintf('\n  %-6s %-26s %-13s %-6s %-7s %-12s %-7s %-14s %-11s %-10s %-10s\n', ...
    'T[N]','tag','verdict','nFlat','sensSt','method','robust', ...
    'red[p n z]','redMinEig','drift','kkt.stat');
for Tn = thrustList(:).'
    if ~isKey(tagOf, Tn)
        fprintf('  %-6g <no tagOf entry -- skipping>\n', Tn);
        continue;
    end
    tag = tagOf(Tn);
    mp  = fullfile(resDir, [tag '.mat']);
    if ~isfile(mp)
        fprintf('  %-6g %-26s MISSING (no cache at %s)\n', Tn, tag, mp);
        continue;
    end
    try
        sosc = verify_sosc_mee(mp);
    catch ME
        fprintf('  %-6g %-26s ERROR (exception: %s)\n', Tn, tag, ME.message);
        T(end+1) = struct('thrustN',Tn,'tag',tag,'verdict','ERROR', ...
            'status','certified-feasibility+sosc-inconclusive','nFlat',NaN, ...
            'sensStable',false,'method','','robust',false, ...
            'red',struct('npos',NaN,'nneg',NaN,'nzero',NaN), ...
            'redMinEig',NaN,'drift',NaN,'stat',NaN); %#ok<AGROW>
        continue;
    end

    save(fullfile(sideDir, ['sosc_' tag '.mat']), 'sosc');

    % Guard the ERROR early-return shape (sosc.red=[] and sosc.kkt=[] when the
    % warm re-solve fails): normalize to NaN placeholders so the row still
    % prints and is appended (an ERROR is a reportable outcome, not a crash).
    red = sosc.red;
    if ~isstruct(red) || isempty(red)
        red = struct('npos',NaN,'nneg',NaN,'nzero',NaN);
    end
    statVal = NaN;
    if isstruct(sosc.kkt) && isfield(sosc.kkt,'stat'), statVal = sosc.kkt.stat; end

    redStr = sprintf('[%g %g %g]', red.npos, red.nneg, red.nzero);
    fprintf('  %-6g %-26s %-13s %-6g %-7d %-12s %-7d %-14s %-11.2e %-10.2e %-10.2e\n', ...
        Tn, tag, sosc.verdict, sosc.nFlat, sosc.sensStable, sosc.method, ...
        sosc.robust, redStr, sosc.redMinEig, sosc.drift, statVal);

    T(end+1) = struct('thrustN',Tn,'tag',tag,'verdict',sosc.verdict, ...
        'status',sosc.status,'nFlat',sosc.nFlat,'sensStable',sosc.sensStable, ...
        'method',sosc.method,'robust',sosc.robust,'red',red, ...
        'redMinEig',sosc.redMinEig,'drift',sosc.drift, ...
        'stat',statVal); %#ok<AGROW>
end
fprintf('\n');
end
