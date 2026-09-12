function P = repolish_endpoints(opts)
%% Purpose:
%
%   RE-POLISH stored rib entries whose costates no longer fly to the endpoint
%   the CURRENT phase rule gives, and re-certify them through the gate stack.
%
%   Why this exists (FINDINGS 50, 2026-09-12). A catalog entry is keyed by
%   its PHASES, not by its endpoint states, so it is only as reproducible as
%   the rule that turns a phase into a state. When `phase_state` replaced a
%   silently-falling-back ordinary spline with the periodic cubic, ONE
%   departure endpoint -- sD = 11/12, the phase nearest the seam on the
%   coarse 105-sample DRO table -- moved by 0.75 mm. Every other grid phase
%   agreed to machine zero. Over a 26.4-day low-thrust arc that 0.75 mm grew
%   into a 2.1 km flown miss, and the foreign witness moved the costates by
%   7.4e-6 to absorb it, which failed the audit.
%
%   The entries are not wrong: each still solves the problem it was
%   certified for. They solve a problem 0.75 mm away from the one the
%   current code poses. The fix is therefore not to re-walk the ribs (hours)
%   but to re-polish each affected entry from its own stored solution
%   against the current endpoints -- a correction of order 1e-5 in the
%   costates, a few Newton steps -- and re-certify it.
%
%   SELF-SELECTING. Nothing is hard-coded to a column. Every stored point is
%   flown from the current endpoints; only those missing by more than
%   opts.tolKm are touched. A run that finds nothing changes nothing.
%
%  ASSUMPTIONS / NOTES:
%
% • Rib files are rewritten IN PLACE after a timestamped backup. The point's
%   provenance fields are preserved and `.repolished` records what happened.
% • A point that fails to re-certify is left AS IT WAS and reported; a
%   half-updated library is worse than a stale one.
% • The sheet (spine) entries sit at sD = 0, a knot of both interpolants, so
%   they are unaffected by construction -- but they are checked anyway.
%
%% Inputs:
%
%  opts                     struct (optional)
%   .resDir ['results'] where the rib files live, .tolKm [0.5] flown miss
%   above which a point is re-polished, .dryRun [false] measure and report
%   without writing, .only [] restrict to these arrival columns, .K [24]
%   segments for the polish, .pool [capped_pool()]
%
%% Outputs:
%
%  P                        struct                  .rows (struct array:
%                                                   .file .j .sD .sA .missBefore
%                                                   .missAfter .dz .ok .reason)
%                                                   .nChecked .nRepolished
%                                                   .nFailed .files (touched)
%
%% Revision History:
%  M. Casey                                                   (c) 09/12/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
resDir = d('resDir', fullfile(here, 'results'));
tolKm  = d('tolKm', 0.5);   dryRun = d('dryRun', false);   K = d('K', 24);
only   = d('only', []);
if isfield(opts, 'pool'), pool = opts.pool; else, pool = capped_pool(); end

lStar = 389703.264829278;
[B, ~] = arclength_arrival('setup');

rows = struct('file', {}, 'j', {}, 'sD', {}, 'sA', {}, 'missBefore', {}, ...
              'missAfter', {}, 'dz', {}, 'ok', {}, 'reason', {});
files = {};
nChecked = 0;  nRe = 0;  nFail = 0;

dd = dir(fullfile(resDir, 'arrival_rib*.mat'));
fprintf('REPOLISH: %d rib file(s) in %s, tolerance %.3f km%s\n', numel(dd), resDir, tolKm, ...
        tern(dryRun, '  (DRY RUN)', ''));

for m = 1:numel(dd)
    f = fullfile(resDir, dd(m).name);
    L = load(f);  R = L.R;  touched = false;
    for q = 1:numel(R)
        if ~isempty(only) && ~ismember(R(q).j, only), continue, end
        pts = R(q).pts;
        % a struct-array element can only be assigned a struct with the SAME
        % fields, so the provenance field is added to the whole array first
        if ~isempty(pts) && ~isfield(pts, 'repolished'), [pts.repolished] = deal([]); end
        for k = 1:numel(pts)
            p = pts(k);
            if ~(isfield(p, 'ok') && ~isempty(p.ok) && p.ok), continue, end
            nChecked = nChecked + 1;
            rv0 = B.stateD(p.sD);  rvf = B.stateA(p.sA);
            miss0 = flownMiss(p.z, rv0, rvf, B, lStar);
            if ~(miss0 >= tolKm), continue, end

            r = struct('file', dd(m).name, 'j', R(q).j, 'sD', p.sD, 'sA', p.sA, ...
                       'missBefore', miss0, 'missAfter', NaN, 'dz', NaN, ...
                       'ok', false, 'reason', '');
            fprintf('  (%s j=%d) sD %.4f sA %.4f: misses %.4f km -> re-polishing\n', ...
                    dd(m).name, R(q).j, p.sD, p.sA, miss0);
            if dryRun
                r.reason = 'dry run';  rows(end+1) = r; %#ok<AGROW>
                continue
            end
            % re-polish from the entry's OWN solution: its junction states if
            % it carries them, else rebuilt from z8
            if isfield(p, 'Y') && ~isempty(p.Y) && size(p.Y, 1) == 14
                seed = struct('tf', p.z(8), 'tGrid', linspace(0, p.z(8), size(p.Y,2)), 'Y', p.Y);
                if size(seed.Y, 2) == K, seed.Y = [seed.Y, seed.Y(:,end)];  seed.tGrid = linspace(0, p.z(8), K+1); end
            else
                seed = seed_from_z8(p.z(:), rv0(1:6), K, B.Tnd, B.cnd, B.mu);
            end
            C = certify_root(seed, rv0(1:6), rvf(1:6), B, ...
                             struct('sA', p.sA, 'sD', p.sD, 'pool', pool));
            if ~C.ok
                r.reason = C.reason;  nFail = nFail + 1;
                fprintf('     REFUSED, entry left unchanged: %s\n', C.reason);
                rows(end+1) = r; %#ok<AGROW>
                continue
            end
            r.dz = norm(C.z(:) - p.z(:));
            r.missAfter = flownMiss(C.z, rv0, rvf, B, lStar);
            r.ok = true;
            % keep every provenance field the point already had; replace the
            % solution and record what was done
            newp = pts(k);                    % same fields and order as the array
            fset = {'z','Y','tfDays','dvKms','mfKg','flyKm','flyVms','dz','conj','g','normR','reason','ok'};
            for ff = fset
                if isfield(C, ff{1}) && isfield(newp, ff{1}), newp.(ff{1}) = C.(ff{1}); end
            end
            newp.ok = true;
            newp.repolished = struct('when', char(datetime('now')), 'missBefore', miss0, ...
                                     'missAfter', r.missAfter, 'dz', r.dz, ...
                                     'why', 'endpoint rule changed (FINDINGS 50)');
            pts(k) = newp;
            touched = true;  nRe = nRe + 1;
            fprintf('     re-certified: |dz| %.3e, miss %.4f -> %.4f km, t_f %.6f d\n', ...
                    r.dz, miss0, r.missAfter, C.tfDays);
            rows(end+1) = r; %#ok<AGROW>
        end
        R(q).pts = pts;
    end
    if touched && ~dryRun
        bak = [f '.bak_endpoint_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))];
        copyfile(f, bak);
        save(f, 'R', '-v7.3');
        files{end+1} = dd(m).name; %#ok<AGROW>
        fprintf('  %s rewritten (backup %s)\n', dd(m).name, bak);
    end
end

P = struct('rows', rows, 'nChecked', nChecked, 'nRepolished', nRe, ...
           'nFailed', nFail, 'files', {files});
fprintf('REPOLISH DONE: %d point(s) checked, %d re-polished, %d refused, %d file(s) rewritten\n', ...
        nChecked, nRe, nFail, numel(files));
end

% ------------------------------------------------------------------------
function km = flownMiss(z8, rv0, rvf, B, lStar)
% FLOWNMISS  Position miss of stored costates flown from the CURRENT
% endpoints.  INPUTS: z8; rv0; rvf; B; lStar.  OUTPUTS: km.
z8 = z8(:);
try
    [~, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6); 1; z8(1:7)], B.Tnd, B.cnd, B.mu);
    km = norm(Y(end,1:3).' - rvf(1:3))*lStar;
catch
    km = Inf;
end
end

% ------------------------------------------------------------------------
function v = tern(c, x, y)
% TERN  Inline conditional.  INPUTS: c; x; y.  OUTPUTS: v.
if c, v = x; else, v = y; end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
