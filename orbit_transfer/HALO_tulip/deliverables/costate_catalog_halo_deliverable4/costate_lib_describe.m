function D = costate_lib_describe(libIn)
%% Purpose:
%
%   Describes a costate library: point it at any of the shipped artifacts --
%   the single-pair library (v1), the thrust-axis library (v2), or the
%   multi-orbit CATALOG -- and it prints the relevant facts a user needs
%   before querying it:
%
%     - the departure-phase grid, in nondimensional time AND days
%     - the arrival-phase grid, in nondimensional time AND days
%     - the thrust levels, Isp and spacecraft mass
%     - the DRO period(s) and tulip petal count(s) covered
%     - coverage statistics: percent of phase pairs solved (per sheet and
%       overall), and which sheet is BEST and which is WORST
%
%   The format is detected automatically; a summary struct is returned for
%   programmatic use.
%
%% Inputs:
%
%  libIn                    struct or char          A loaded library/catalog
%                                                   struct, or the path of a
%                                                   .mat holding one
%
%% Outputs:
%
%  D                        struct                  .format ('v1'|'v2'|
%                                                   'catalog'), .nEntries,
%                                                   .coveragePct (overall),
%                                                   .sheets (per-sheet
%                                                   coverage, catalog only),
%                                                   .best / .worst sheet
%                                                   labels
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: describe the catalog sitting beside this function's caller:
       D = costate_lib_describe('costate_catalog_dro_tulip.mat');
     return;
end

%% Accept a filename or a struct:
if ischar(libIn) || isstring(libIn)
    S = load(char(libIn));
    fn = fieldnames(S);
    lib = S.(fn{1});
else
    lib = libIn;
end

%% Detect the format:
if isfield(lib, 'sheets'),        fmt = 'catalog';
elseif ndims(lib.grid.has_solution) == 3 && isfield(lib.grid,'thrust_N')
                                  fmt = 'v2';
else,                             fmt = 'v1';
end
tStar = pickfield(lib, {'constants.tStar_s'}, 382981.289129055);
d2day = tStar/86400;

fprintf('\n==================== COSTATE LIBRARY ====================\n');
fprintf('%s  (%s format)\n', lib.name, fmt);
fprintf('created %s;  %d entries\n', lib.created, lib.n_entries);

switch fmt
%% ------------------------------------------------------------------ v1/v2
case {'v1','v2'}
    Pd = lib.departure_params.period_days;
    Pa = lib.arrival_params.period_days;
    sD = lib.grid.departure_phase_frac;
    sA = lib.grid.arrival_phase_frac;
    fprintf('\nDRO period: %.4f ND = %.3f days (tau = %.3f)\n', ...
            Pd/d2day, Pd, lib.departure_params.tau);
    fprintf('tulip: Np = %d, period %.4f ND = %.3f days\n', ...
            lib.arrival_params.Np, Pa/d2day, Pa);
    fprintf('\ndeparture phases (%d): ND  [%s]\n', numel(sD), ...
            sprintf('%.3f ', sD*Pd/d2day));
    fprintf('                      days [%s]\n', sprintf('%.2f ', sD*Pd));
    fprintf('arrival phases (%d):   ND  [%s]\n', numel(sA), ...
            sprintf('%.3f ', sA*Pa/d2day));
    fprintf('                      days [%s]\n', sprintf('%.2f ', sA*Pa));
    if strcmp(fmt,'v2')
        rungs = lib.grid.thrust_N;
        fprintf('\nthrust levels (N): [%s] at Isp %d s, m0 %d kg\n', ...
                num2str(rungs), lib.thruster.isp_s, lib.thruster.m0_kg);
        anyPair = any(lib.grid.has_solution, 3);
        fprintf('\ncoverage: %d/%d phase pairs solved at >=1 thrust (%.0f%%)\n', ...
                nnz(anyPair), numel(anyPair), 100*nnz(anyPair)/numel(anyPair));
        fprintf('per rung:');
        for kr = 1:numel(rungs)
            fprintf('  %.3g N: %d', rungs(kr), nnz(lib.grid.has_solution(:,:,kr)));
        end
        fprintf('\n');
        cov = nnz(anyPair)/numel(anyPair);
    else
        fprintf('\nthrust: %.3g N at Isp %d s, m0 %d kg\n', ...
                lib.thruster.Tmax_N, lib.thruster.isp_s, lib.thruster.m0_kg);
        cov = nnz(lib.grid.has_solution)/numel(lib.grid.has_solution);
        fprintf('coverage: %d/%d phase pairs solved (%.0f%%)\n', ...
                nnz(lib.grid.has_solution), numel(lib.grid.has_solution), 100*cov);
    end
    D = struct('format',fmt, 'nEntries',lib.n_entries, 'coveragePct',100*cov);

%% ----------------------------------------------------------------- catalog
case 'catalog'
    sh = lib.sheets;
    taus = unique([sh.tauDRO]);
    nps  = unique([sh.Np]);
    fprintf('\nDRO periods (tau = period): ND  [%s]\n', sprintf('%.2f ', taus));
    fprintf('                            days [%s]\n', sprintf('%.2f ', taus*d2day));
    fprintf('tulip petal counts: [%s]  (pm = %d branch; +1 is its z-mirror)\n', ...
            num2str(nps), sh(1).pm);
    fprintf('tulip periods:      ND  [%s]\n', ...
            sprintf('%.3f ', unique([sh.period_tulip_nd])));
    fprintf('                    days [%s]\n', ...
            sprintf('%.2f ', unique([sh.period_tulip_nd])*d2day));
    fprintf('thrust levels (N): [%s] at Isp %d s, m0 %d kg\n', ...
            num2str(lib.rungs_N), lib.thruster.isp_s, lib.thruster.m0_kg);
    sh1 = sh(1);
    fprintf(['\nphasing grid per sheet: %d departure x %d arrival phases\n', ...
             '  (fractions of each orbit period; e.g. sheet 1 departures ND [%s],\n', ...
             '   days [%s])\n'], ...
            numel(sh1.sD_frac), numel(sh1.sA_frac), ...
            sprintf('%.3f ', sh1.sD_frac*sh1.tauDRO), ...
            sprintf('%.2f ', sh1.sD_frac*sh1.tauDRO*d2day));

    %% per-sheet statistics -> best and worst
    nS = numel(sh);
    covPct = zeros(nS,1);  nEnt = zeros(nS,1);  lab = cell(nS,1);
    fprintf('\nsheet                pairs solved   entries   full ladders\n');
    for k = 1:nS
        anyP = any(sh(k).has_solution, 3);
        covPct(k) = 100*nnz(anyP)/numel(anyP);
        nEnt(k) = size(sh(k).z8, 2);
        lab{k} = sprintf('tau=%.2f Np=%2d', sh(k).tauDRO, sh(k).Np);
        fprintf('  %-16s   %2d/%2d (%3.0f%%)   %5d      %d\n', lab{k}, ...
                nnz(anyP), numel(anyP), covPct(k), nEnt(k), ...
                nnz(sum(sh(k).has_solution,3) == numel(lib.rungs_N)));
    end
    % rank by coverage, entries as tiebreak
    score = covPct + nEnt/max(nEnt)/100;
    [~, kBest]  = max(score);
    [~, kWorst] = min(score);
    fprintf('\noverall: %.0f%% of phase pairs solved;  %d entries total\n', ...
            mean(covPct), lib.n_entries);
    fprintf('BEST  sheet: %s  (%.0f%% coverage, %d entries)\n', ...
            lab{kBest}, covPct(kBest), nEnt(kBest));
    fprintf('WORST sheet: %s  (%.0f%% coverage, %d entries)\n', ...
            lab{kWorst}, covPct(kWorst), nEnt(kWorst));
    D = struct('format',fmt, 'nEntries',lib.n_entries, ...
               'coveragePct',mean(covPct), ...
               'sheets',struct('label',{lab},'coveragePct',covPct,'entries',nEnt), ...
               'best',lab{kBest}, 'worst',lab{kWorst});
end
fprintf('=========================================================\n\n');
end

% ------------------------------------------------------------------------
function v = pickfield(s, paths, vDef)
% PICKFIELD  First existing nested field among paths, else the default.
% INPUTS: s struct; paths cell of 'a.b' strings; vDef default. OUTPUTS: v.
v = vDef;
for k = 1:numel(paths)
    parts = strsplit(paths{k}, '.');
    t = s;  ok = true;
    for q = 1:numel(parts)
        if isfield(t, parts{q}), t = t.(parts{q}); else, ok = false; break, end
    end
    if ok, v = t; return, end
end
end
