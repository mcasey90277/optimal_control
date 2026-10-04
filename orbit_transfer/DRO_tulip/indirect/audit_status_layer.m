function A_ = audit_status_layer(catMat, opts)
%% Purpose:
%
%   The FAIL-CLOSED audit of a catalog's optimality-status layer. Every
%   primary must carry status 4 (sufficient) and junctions that, flown
%   segment by segment (entry_thrust_program), join up and reach the
%   arrival orbit. Every alternative is re-certified from its OWN junctions:
%   the stored status must be reproduced, and the root must not move -- a
%   re-certification that lands on another root is BAD, never relabelled
%   (spec 7, Review Focus 5). Anything not shown to hold is a problem.
%
%  ASSUMPTIONS / NOTES:
%
% • One sheet per catalog. The endpoint closures B.stateD / B.stateA work
%   at any phase, so one arclength_arrival setup serves every entry.
% • Junctions sit on linspace(0, t_f, K+1); column k is the augmented state
%   at the start of segment k. The alternatives' seed repeats the last
%   column as the K+1 node and forces node 1 to [rv0; 1; z8(1:7)].
% • .certifier is a TEST SEAM only; production uses certify_root with a
%   pool ('allowUnfenced' is set only when no pool is passed).
% • A certifier that throws is a BAD row, not a skipped one.
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
%                                                   .tolMove [1e-6] (same_root
%                                                   carries its own 1e-6)
%                                                   .junctionKm [1] max defect
%                                                   .flyKm [1] arrival miss
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
junctionKm    = fieldd(opts, 'junctionKm', 1);
flyKm         = fieldd(opts, 'flyKm', 1);
certifier     = fieldd(opts, 'certifier', @certify_root);
if isfield(opts, 'pool'), pool = opts.pool; else, pool = capped_pool(); end

L = load(catMat);  c = L.(char(fieldnames(L)));
assert(isscalar(c.sheets), 'audit_status_layer:sheets', 'one sheet expected');
sh = c.sheets(1);
problems = {};
primRows = struct('k', {}, 'ok', {}, 'why', {});
altRows  = struct('k', {}, 'ok', {}, 'why', {}, 'statusNow', {}, 'moved', {});
B = [];

% ---- primaries: status 4 and junctions that fly ----------------------------
if ~skipPrimaries
    assert(isfield(sh, 'status') && isfield(sh, 'junctions'), 'audit_status_layer:noLayer', ...
           'the catalog has no status layer (status, junctions): nothing to audit');
    [B, phys] = setupAt(c, sh.sD_frac(1));
    for iA = 1:numel(sh.sA_frac)
        for iD = 1:numel(sh.sD_frac)
            if ~sh.has_solution(iD, iA, 1), continue, end
            k = sh.entry_index(iD, iA, 1);  why = '';
            if sh.status(iD, iA, 1) ~= 4, why = sprintf('primary status %d, not 4', sh.status(iD, iA, 1)); end
            Jn = [];  if isempty(why) && k <= numel(sh.junctions), Jn = sh.junctions{k}; end
            if isempty(why) && isempty(Jn), why = 'no junctions'; end
            if isempty(why)
                try
                    rv0 = B.stateD(sh.sD_frac(iD));  rvf = B.stateA(sh.sA_frac(iA));
                    Pp = entry_thrust_program(sh.z8(:, k), Jn, rv0(1:6), rvf(1:6), phys);
                    if max([Pp.defectKm 0]) > junctionKm, why = sprintf('junction defect %.2f km', max(Pp.defectKm));
                    elseif ~(Pp.flyKm <= flyKm), why = sprintf('junction flight misses by %.2f km', Pp.flyKm);
                    elseif ~(Pp.startErr < 1e-9), why = sprintf('first junction differs from z8 (%.1e)', Pp.startErr); end
                catch ME
                    why = ['junction flight threw: ' ME.message];
                end
            end
            primRows(end+1) = struct('k', k, 'ok', isempty(why), 'why', why);
            if ~isempty(why), problems{end+1} = sprintf('primary (%d,%d): %s', iD, iA, why); end
        end
    end
end

% ---- alternatives: re-certify from their own junctions ---------------------
Alt = struct([]);  if isfield(c, 'alternatives'), Alt = c.alternatives; end
idx = 1:numel(Alt);  if ~isempty(idxAlt), idx = idxAlt(:).'; end
if ~isempty(idx) && isempty(B), [B, ~] = setupAt(c, sh.sD_frac(1)); end
for k = idx
    a = Alt(k);  K = size(a.junctions, 2);
    try
        rv0 = B.stateD(a.sD);  rvf = B.stateA(a.sA);
        seed = struct('tf', a.z8(8), 'tGrid', linspace(0, a.z8(8), K + 1), 'Y', [a.junctions, a.junctions(:, end)]);
        seed.Y(1:7, 1) = [rv0(1:6); 1];  seed.Y(8:14, 1) = a.z8(1:7);
        C = certifier(seed, rv0(1:6), rvf(1:6), B, struct('sD', a.sD, 'sA', a.sA, 'pool', pool, 'allowUnfenced', isempty(pool)));
        moved = ~same_root(C.z, a.z8);
        if moved, why = 'the re-certification moved to another root';
        elseif C.status ~= a.status, why = sprintf('status not reproduced: stored %d, now %d (%s)', a.status, C.status, C.reason);
        else, why = ''; end
        sNow = C.status;
    catch ME
        moved = false;  sNow = NaN;  why = ['re-certification threw: ' ME.message];
    end
    altRows(end+1) = struct('k', k, 'ok', isempty(why), 'why', why, 'statusNow', sNow, 'moved', moved);
    if ~isempty(why), problems{end+1} = sprintf('alternative %d (%.4f, %.4f): %s', k, a.sD, a.sA, why); end
    if ~isempty(outFile), save(outFile, 'primRows', 'altRows', 'problems'); end    % resume-friendly
end

% ---- verdict -----------------------------------------------------------------
A_ = struct('primRows', primRows, 'altRows', altRows, ...
            'nOk', sum([primRows.ok]) + sum([altRows.ok]), 'nBad', numel(problems), ...
            'contentKey', catalog_content_key(c), 'altContentKey', alternatives_content_key(c));
A_.problems = problems;
if ~isempty(outFile), audit = A_; save(outFile, 'audit', 'primRows', 'altRows', 'problems'); end
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
