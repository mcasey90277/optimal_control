function A_ = audit_status_layer(catMat, opts)
%% Purpose:
%
%   The FAIL-CLOSED audit of a catalog's optimality-status layer. Every
%   primary must carry status 4 (sufficient) and junctions that, flown
%   segment by segment (entry_thrust_program), join up in position AND
%   velocity and reach the arrival orbit. Every alternative must first fly
%   from its own junctions (junctions join up, arrival within the library
%   floor), then be re-certified from those junctions: the stored status
%   must be reproduced, and the root must not move -- a re-certification
%   that lands on another root is BAD, never relabelled (spec 7, Review
%   Focus 5). Anything not shown to hold is a problem.
%
%  ASSUMPTIONS / NOTES:
%
% • One sheet per catalog. The endpoint closures B.stateD / B.stateA work
%   at any phase, so one arclength_arrival setup serves every entry. If the
%   setup fails, every row it would have served is BAD with the reason and
%   the audit returns normally.
% • Junctions sit on linspace(0, t_f, K+1); column k is the augmented state
%   at the start of segment k. The alternatives' seed repeats the last
%   column as the K+1 node and forces node 1 to [rv0; 1; z8(1:7)].
% • Why a FULL certifier re-run per alternative: certify_root returns at
%   its first failing gate, so reproducing the stored status from the
%   row's own junctions is itself the targeted check -- a status that no
%   longer holds stops at the gate that broke and reports it.
% • "Same root" is same_root's rule for the whole library (1e-6); the audit
%   has no tolerance of its own to loosen it.
% • The fence is never opened implicitly: with no pool, alternatives are
%   BAD ("no pool: refusing to re-certify unfenced") unless the caller sets
%   .allowUnfenced (the seam test does). .certifier is a TEST SEAM only.
% • A throw anywhere in a row (setup, flight, certifier) is a BAD row, and
%   every tolerance test is written so that NaN fails.
% • .out: the partial save after each alternative is CRASH SALVAGE (rows so
%   far + both content keys), not a resume.
%
%% Inputs:
%
%  catMat                   char                    catalog .mat (one variable)
%
%  opts                     struct (optional)       .idxAlt [] alternatives to
%                                                   audit (default all)
%                                                   .skipPrimaries [false]
%                                                   .out '' save path
%                                                   .pool [capped_pool()]
%                                                   .allowUnfenced [false]
%                                                   .junctionKm [1] max junction
%                                                   position defect (km)
%                                                   .junctionVms [1] max junction
%                                                   velocity defect (m/s)
%                                                   .flyKm [1] .flyVms [0.1]
%                                                   primary arrival miss
%                                                   .floorKm [100] .floorVms [10]
%                                                   alternative arrival miss
%                                                   (the library floor)
%                                                   .certifier [] @(seed, rv0,
%                                                   rvf, B, opts) -> C (.z
%                                                   .status .reason)
%
%% Outputs:
%
%  A_                       struct                  .primRows (.k .ok .why)
%                                                   .altRows (.k .ok .why
%                                                   .statusNow .moved) .nOk
%                                                   .nBad .problems {1 x n}
%                                                   .contentKey .altContentKey
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
idxAlt        = fieldd(opts, 'idxAlt', []);
skipPrimaries = fieldd(opts, 'skipPrimaries', false);
outFile       = fieldd(opts, 'out', '');
allowUnfenced = fieldd(opts, 'allowUnfenced', false);
tolPrim       = struct('defKm', fieldd(opts, 'junctionKm', 1), 'defVms', fieldd(opts, 'junctionVms', 1), ...
                       'flyKm', fieldd(opts, 'flyKm', 1), 'flyVms', fieldd(opts, 'flyVms', 0.1));
tolAlt        = tolPrim;                                     % alternatives: the floor, not the primary gate
tolAlt.flyKm  = fieldd(opts, 'floorKm', 100);
tolAlt.flyVms = fieldd(opts, 'floorVms', 10);
certifier     = fieldd(opts, 'certifier', @certify_root);
if isfield(opts, 'pool'), pool = opts.pool; else, pool = capped_pool(); end

L = load(catMat);  c = L.(char(fieldnames(L)));
assert(isscalar(c.sheets), 'audit_status_layer:sheets', 'one sheet expected');
sh = c.sheets(1);
contentKey = catalog_content_key(c);  altContentKey = alternatives_content_key(c);
problems = {};
primRows = struct('k', {}, 'ok', {}, 'why', {});
altRows  = struct('k', {}, 'ok', {}, 'why', {}, 'statusNow', {}, 'moved', {});
Alt = struct([]);  if isfield(c, 'alternatives'), Alt = c.alternatives; end
idx = 1:numel(Alt);  if ~isempty(idxAlt), idx = idxAlt(:).'; end

% ---- one setup for every row ---------------------------------------------------
setupWhy = '';  B = [];  phys = [];
if ~skipPrimaries || ~isempty(idx)
    try
        [B, phys] = setupAt(c, sh.sD_frac(1));
    catch ME
        setupWhy = ['setup failed: ' ME.message];
    end
end

% ---- primaries: status 4 and junctions that fly ----------------------------
if ~skipPrimaries
    assert(isfield(sh, 'status') && isfield(sh, 'junctions'), 'audit_status_layer:noLayer', ...
           'the catalog has no status layer (status, junctions): nothing to audit');
    for iA = 1:numel(sh.sA_frac)
        for iD = 1:numel(sh.sD_frac)
            if ~sh.has_solution(iD, iA, 1), continue, end
            k = sh.entry_index(iD, iA, 1);  why = setupWhy;
            if isempty(why) && sh.status(iD, iA, 1) ~= 4, why = sprintf('primary status %d, not 4', sh.status(iD, iA, 1)); end
            Jn = [];  if isempty(why) && k <= numel(sh.junctions), Jn = sh.junctions{k}; end
            if isempty(why) && isempty(Jn), why = 'no junctions'; end
            if isempty(why)
                why = flightWhy(sh.z8(:, k), Jn, B, sh.sD_frac(iD), sh.sA_frac(iA), phys, tolPrim);
            end
            primRows(end+1) = struct('k', k, 'ok', isempty(why), 'why', why);
            if ~isempty(why), problems{end+1} = sprintf('primary (%d,%d): %s', iD, iA, why); end
        end
    end
end

% ---- alternatives: fly, then re-certify from their own junctions ----------
for k = idx
    a = Alt(k);  moved = false;  sNow = NaN;  why = setupWhy;
    if isempty(why) && isempty(a.junctions), why = 'no junctions'; end
    if isempty(why), why = flightWhy(a.z8, a.junctions, B, a.sD, a.sA, phys, tolAlt); end
    if isempty(why) && isempty(pool) && ~allowUnfenced, why = 'no pool: refusing to re-certify unfenced'; end
    if isempty(why)
        try
            K = size(a.junctions, 2);
            rv0 = B.stateD(a.sD);  rvf = B.stateA(a.sA);
            seed = struct('tf', a.z8(8), 'tGrid', linspace(0, a.z8(8), K + 1), 'Y', [a.junctions, a.junctions(:, end)]);
            seed.Y(1:7, 1) = [rv0(1:6); 1];  seed.Y(8:14, 1) = a.z8(1:7);
            C = certifier(seed, rv0(1:6), rvf(1:6), B, struct('sD', a.sD, 'sA', a.sA, 'pool', pool, ...
                                                              'allowUnfenced', allowUnfenced));
            moved = ~same_root(C.z, a.z8);
            if moved, why = 'the re-certification moved to another root';
            elseif ~(C.status == a.status), why = sprintf('status not reproduced: stored %d, now %d (%s)', a.status, C.status, C.reason); end
            sNow = C.status;
        catch ME
            moved = false;  sNow = NaN;  why = ['re-certification threw: ' ME.message];
        end
    end
    altRows(end+1) = struct('k', k, 'ok', isempty(why), 'why', why, 'statusNow', sNow, 'moved', moved);
    if ~isempty(why), problems{end+1} = sprintf('alternative %d (%.4f, %.4f): %s', k, a.sD, a.sA, why); end
    if ~isempty(outFile), save(outFile, 'primRows', 'altRows', 'problems', 'contentKey', 'altContentKey'); end   % crash salvage
end

% ---- verdict -----------------------------------------------------------------
A_ = struct('primRows', primRows, 'altRows', altRows, ...
            'nOk', sum([primRows.ok]) + sum([altRows.ok]), 'nBad', numel(problems), ...
            'contentKey', contentKey, 'altContentKey', altContentKey);
A_.problems = problems;
if ~isempty(outFile), audit = A_; save(outFile, 'audit', 'primRows', 'altRows', 'problems', 'contentKey', 'altContentKey'); end
end

function why = flightWhy(z8, Jn, B, sD, sA, phys, tol)
% FLIGHTWHY  Fly an entry from its junctions; '' when every check holds.
% INPUTS:  z8 [8x1]; Jn [14xK] junctions; B [struct] setup; sD, sA [scalar]
%          phases; phys [struct]; tol [struct] .defKm .defVms .flyKm .flyVms.
% OUTPUTS: why [char] first failed check, '' when none (NaN fails).
why = '';
try
    rv0 = B.stateD(sD);  rvf = B.stateA(sA);
    P = entry_thrust_program(z8, Jn, rv0(1:6), rvf(1:6), phys);
    if ~all(P.defectKm <= tol.defKm), why = sprintf('junction defect %.2f km', max(P.defectKm));
    elseif ~all(P.defectVms <= tol.defVms), why = sprintf('junction defect %.3f m/s', max(P.defectVms));
    elseif ~(P.flyKm <= tol.flyKm), why = sprintf('junction flight misses by %.2f km', P.flyKm);
    elseif ~(P.flyVms <= tol.flyVms), why = sprintf('junction flight misses by %.4f m/s', P.flyVms);
    elseif ~(P.startErr < 1e-9), why = sprintf('first junction differs from z8 (%.1e)', P.startErr); end
catch ME
    why = ['junction flight threw: ' ME.message];
end
end

function [B, phys] = setupAt(c, sD)
% SETUPAT  Endpoint closures and propulsion constants for the catalog's problem.
% INPUTS:  c [struct] catalog; sD [scalar] departure phase.
% OUTPUTS: B [struct] arclength_arrival setup; phys [struct] .Tnd .cnd .mu
%          .lStar .tStar .m0kg for entry_thrust_program.
[B, ~] = arclength_arrival('setup', catalog_setup_request(c, sD));
phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'mu', B.mu, 'lStar', c.constants.lStar_km, ...
              'tStar', c.constants.tStar_s, 'm0kg', c.thruster.m0_kg);
end

function v = fieldd(s, f, d_)
% FIELDD  Field of s, or the default when absent or empty.
% INPUTS:  s [struct]; f [char]; d_ default.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
