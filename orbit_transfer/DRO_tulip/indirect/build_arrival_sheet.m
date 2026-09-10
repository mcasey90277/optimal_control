function S = build_arrival_sheet(opts)
%% Purpose:
%
%   Front door for the ARRIVAL-PHASE sheet at one operating point: gather
%   every arclength arc saved in results/ (arrival_arc_*.mat), certify all
%   grid crossings with certify_crossing, assemble the sheet with
%   sheet_from_arcs, save it, and print the table (grid phase, minimum
%   certified t_f, dV, fuel, number of candidates / certified, first
%   failure reason when nothing certified).
%
%% Inputs:
%
%  opts                     struct (optional)
%   .pattern ['arrival_arc_*.mat'] .out ['results/arrival_sheet_70mN.mat']
%   .sA0 [0.0754] .nA [12] .copts (certify_crossing options)
%   plus arclength_arrival setup options (.thrustN .ispS .sD ...)
%
%% Outputs:
%
%  S                        struct                  sheet_from_arcs output +
%                                                   .arcs (file list) .B
%                                                   .anc .opts .built
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));
pat = d('pattern', 'arrival_arc_*.mat');
out = d('out', fullfile(here, 'results', 'arrival_sheet_70mN.mat'));
lStar = 389703.264829278;  tStar = 382981.289129055;

% the fence's pool, created AFTER startup has set the path (workers inherit
% the client path at pool creation)
pool = capped_pool();
[B, anc] = arclength_arrival('setup', opts);
files = dir(fullfile(here, 'results', pat));
assert(~isempty(files), 'no arcs match %s', pat);
arcs = cell(1, numel(files));
for k = 1:numel(files)
    L = load(fullfile(files(k).folder, files(k).name));
    arcs{k} = L.A;
    fprintf('arc %d: %-32s %4d roots, sA %.4f -> %.4f, %d folds, %d crossings, stop = %s\n', ...
        k, files(k).name, numel(L.A.q), L.A.q(1), L.A.q(end), numel(L.A.folds), numel(L.A.crossings), L.A.stop);
end

% SEED the sheet with the certified library solutions at this departure
% phase: they are the arcs' start points, so they are not crossings of any
% arc and a crossings-only sheet reports NaN at exactly the phases whose
% solutions launched it.
sD0 = anc.sD;  nA = d('nA', 12);  sA0 = d('sA0', 0.0754);
policy = d('copts', struct());
policy.pool = pool;
if ~isfield(policy, 'wallSec'), policy.wallSec = 600; end
lib = dro_tulip_library(here);
lib = lib(abs(mod([lib.sD] - sD0 + 0.5, 1) - 0.5) < 1e-8);
seeds = struct([]);
for k = 1:numel(lib)
    g = mod(lib(k).sA - sA0, 1)*nA;
    if abs(g - round(g)) > 1e-6 && abs(g - nA) > 1e-6, continue, end   % off-grid
    K = size(lib(k).Y, 2);
    seed = struct('tf', lib(k).z(8), 'tGrid', linspace(0, lib(k).z(8), K+1), ...
                  'Y', [lib(k).Y, lib(k).Y(:,end)]);
    rv0 = B.stateD(sD0);
    seed.Y(1:7,1) = [rv0(1:6); 1];  seed.Y(8:14,1) = lib(k).z(1:7);
    % ONE certification policy for both routes. Seeds used to be certified
    % with a fresh hardcoded struct while opts.copts applied only to
    % crossings, so a stricter requested gate silently did not reach the
    % seeds that populate the same sheet. (Astra chain review 2026-09-10.)
    C = certify_root(seed, rv0, B.stateA(lib(k).sA), B, ...
                     setfield(setfield(policy, 'sA', lib(k).sA), 'sD', sD0)); %#ok<SFLD>
    fprintf('library seed (%.4f, %.4f) [%s]: %s\n', sD0, lib(k).sA, lib(k).src, C.reason);
    if isempty(seeds), seeds = C; else, seeds(end+1) = C; end %#ok<AGROW>
end

S = sheet_from_arcs(arcs, struct('sA0', sA0, 'nA', nA, 'seeds', seeds, ...
                                 'B', B, 'anc', anc, 'copts', policy));
S.arcs = {files.name};  S.B = B;  S.anc = anc;  S.opts = opts;  S.built = datestr(now);
S.problem = B.problem;          % the identity packaging must use
S.policy = rmfield(policy, 'pool');
S.B = rmfield(S.B, {'res', 'dRdq', 'stateA', 'stateD'});
save(out, 'S');

fprintf('\n  j    sA      t_f [d]   dV [km/s]  fuel [kg]  cand  cert   note\n');
for j = 1:numel(S.sA)
    c = S.cand{j};  nc = numel(c);  ncert = 0;  note = '';
    if nc > 0, ncert = nnz([c.ok]); end
    if isfinite(S.TF(j))
        k = find([c.ok] & abs([c.tfDays] - S.TF(j)) < 1e-9, 1);
        fprintf(' %2d  %.4f  %9.4f  %9.4f  %9.2f   %2d    %2d\n', j, S.sA(j), S.TF(j), c(k).dvKms, c(k).mfKg, nc, ncert);
    else
        if nc > 0, note = c(1).reason; else, note = 'no crossing'; end
        fprintf(' %2d  %.4f        ---        ---        ---   %2d    %2d   %s\n', j, S.sA(j), nc, ncert, note);
    end
end
fprintf('saved %s\n', out);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
