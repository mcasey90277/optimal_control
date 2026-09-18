function A = audit_phase_catalog(catMat, opts)
%% Purpose:
%
%   Audit a SHIPPED phase catalog the way a RECIPIENT would: rebuild each
%   entry's endpoints from the catalog's own keys, fly its own z8, and check
%   that the trajectory does what the catalog says it does. Nothing internal
%   to the build is consulted -- if the catalog is self-consistent and
%   truthful, this passes on the file alone.
%
%   Per entry it checks:
%     LABEL   the stored t_f agrees with z8(8); the phases lie on the grid
%     FLIGHT  flying z8 from the labelled departure phase reaches the
%             labelled arrival phase (position and velocity), completes to
%             t_f, and obeys the all-burn mass law
%     WITNESS pumpkyn tfMin from z8 returns z8, and ITS solution also flies
%             to the target
%     SECOND  the conjugate verdict and the three hypothesis gates,
%             RECOMPUTED and compared with what the catalog stores. FAIL
%             CLOSED: a re-polish or gates call that fails or times out, a
%             polished root that left the stored one (.tolMove), a gate
%             that is not a pass (H2, H3, dim S = 1, H6) -- each is a BAD
%             row, never a skipped check
%
%   Written 2026-09-10 after an external review found labelling and export
%   defects in the builder. The builder has been fixed; this answers the
%   separate question of whether the ALREADY SHIPPED file is correct, which
%   no amount of fixing the builder can answer.
%
%% Inputs:
%
%  catMat                   char                    catalog .mat
%  opts                     struct (optional)
%   .gateKm [100] .gateVms [10] gross flight screens, .tolFlyKm [1] the
%   ENDPOINT-REPRODUCTION gate (see below), .tolDz [1e-6] .idx [] entries to audit
%   (default all), .out '' save path, .pool [capped_pool()], .tolMove [1e-6],
%   .hypFloor [1e-5] the certifier's H2/H3 floor on the whole-arc lower bounds,
%   .cappedWrap [] @(realCap) -> a replacement fence (tests/test_audit_fail_closed)
%
%% Outputs:
%
%  A                        struct                  .rows (per entry) .nBad
%                                                   .nOk .problems (cellstr)
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));
gateKm = d('gateKm', 100);  gateVms = d('gateVms', 10);  tolDz = d('tolDz', 1e-6);
% ENDPOINT REPRODUCTION. A stored entry is keyed by its PHASES, so it is only
% as reproducible as the rule that turns a phase into a state. When that rule
% changes, the entry still solves its own problem but no longer solves the
% one the current code poses, and the symptom is a flown miss that grows far
% beyond propagation error. Measured 2026-09-12: replacing a silent ordinary
% spline with the periodic cubic moved ONE departure endpoint by 0.75 mm and
% the flown miss of its 26.4-day arc from 0.064 km to 2.125 km, which
% surfaced downstream as an unexplained witness disagreement. This gate names
% it at the point of occurrence. It sits well above honest propagation error
% (the worst legitimate entry measures 0.29 km) and well below that failure.
tolFlyKm = d('tolFlyKm', 1);
pool = d('pool', capped_pool());
% the fence. A test may WRAP it (opts.cappedWrap = @(realCap) wrappedCap) to
% knock one call out; the real one is handed in because it is local here
cap = @capped;
if isfield(opts, 'cappedWrap') && ~isempty(opts.cappedWrap), cap = opts.cappedWrap(@capped); end
tolMove = d('tolMove', 1e-6);               % how far the re-polish may move a stored root
hypFloor = d('hypFloor', 1e-5);             % H2/H3 floor on the between-sample bounds (as certify_root)

L = load(catMat);  fn = fieldnames(L);  c = L.(fn{1});
assert(numel(c.sheets) == 1, 'this audit handles a single-sheet phase catalog');
s = c.sheets(1);
lStar = c.constants.lStar_km;  tStar = c.constants.tStar_s;  mu = c.constants.muStar;
m0 = c.thruster.m0_kg;  isp = c.thruster.isp_s;
% THE shared propulsion conversion (costate_common/nd_propulsion)
ndp = nd_propulsion(c.rungs_N(1), isp, m0, lStar, tStar);
cnd = ndp.cnd;   Tnd = ndp.Tnd;

% endpoints rebuilt from the catalog's own keys, not from our build state
ob = struct('muStar', mu, 'lStar', lStar, 'tStar', tStar, 'tauDRO', s.tauDRO, ...
            'NpTulip', s.Np, 'tauTulip', s.period_tulip_nd, 'pmTulip', s.pm, ...
            'ispS', isp, 'm0kg', m0);
[tD, rvD, tT, rvT] = ladder_endpoints(ob);
% THE shared endpoint rule (costate_common/phase_state, FINDINGS 44)
stD = phase_state(tD, rvD);   stA = phase_state(tT, rvT);

[iD, iA, iR] = ind2sub(size(s.has_solution), find(s.has_solution));
idx = d('idx', 1:numel(iD));
rows = struct('iD', {}, 'iA', {}, 'sD', {}, 'sA', {}, 'tfStored', {}, 'tfZ8', {}, ...
              'flyKm', {}, 'flyVms', {}, 'dz', {}, 'witKm', {}, 'massErr', {}, ...
              'conjStored', {}, 'conjNow', {}, 'dimSStored', {}, 'dimSNow', {}, ...
              'dzPolish', {}, 'ok', {}, 'problem', {});
problems = {};
for kk = idx(:)'
    r = struct('iD', iD(kk), 'iA', iA(kk), 'sD', s.sD_frac(iD(kk)), 'sA', s.sA_frac(iA(kk)), ...
               'tfStored', s.tf_nd(iD(kk), iA(kk), iR(kk)), 'tfZ8', NaN, 'flyKm', NaN, ...
               'flyVms', NaN, 'dz', NaN, 'witKm', NaN, 'massErr', NaN, ...
               'conjStored', NaN, 'conjNow', NaN, 'dimSStored', NaN, 'dimSNow', NaN, ...
               'dzPolish', NaN, 'ok', false, 'problem', '');
    z8 = s.z8(:, s.entry_index(iD(kk), iA(kk), iR(kk)));
    r.tfZ8 = z8(8);
    if isfield(s, 'conj_pass'),  r.conjStored  = double(s.conj_pass(iD(kk), iA(kk), iR(kk))); end
    if isfield(s, 'gate_dimS'),  r.dimSStored  = s.gate_dimS(iD(kk), iA(kk), iR(kk)); end

    % LABEL
    if abs(r.tfStored - r.tfZ8) > 1e-12*max(r.tfZ8, 1)
        r.problem = sprintf('stored t_f %.12g disagrees with z8(8) %.12g', r.tfStored, r.tfZ8);
        rows(end+1) = r;  problems{end+1} = r.problem; %#ok<AGROW>
        continue
    end

    rv0 = stD(r.sD);  rvf = stA(r.sA);
    % FLIGHT
    [okF, tF, Yf] = cap(pool, 300, @pumpkyn.cr3bp.tfMinProp, 2, z8(8), [rv0(1:6); 1; z8(1:7)], Tnd, cnd, mu);
    if ~okF, r.problem = 'flight timed out';  rows(end+1) = r;  problems{end+1} = r.problem; continue, end %#ok<AGROW>
    r.flyKm  = norm(Yf(end,1:3) - rvf(1:3)')*lStar;
    r.flyVms = norm(Yf(end,4:6) - rvf(4:6)')*lStar/tStar*1000;
    mExp = 1 - (Tnd/cnd)*z8(8);
    r.massErr = abs(Yf(end,7) - mExp);
    bad = {};
    if abs(tF(end) - z8(8)) > 1e-8*max(z8(8),1), bad{end+1} = 'flight did not reach t_f'; end
    if ~(r.flyKm  < gateKm),  bad{end+1} = sprintf('arrival %.1f km', r.flyKm); end
    if r.flyKm >= tolFlyKm && r.flyKm < gateKm
        bad{end+1} = sprintf(['ENDPOINT REPRODUCTION: the stored costates fly %.3f km from the ' ...
                              'endpoint the CURRENT phase rule gives (limit %g km) -- the entry ' ...
                              'was certified against a different endpoint; re-polish it'], ...
                             r.flyKm, tolFlyKm);
    end
    if ~(r.flyVms < gateVms), bad{end+1} = sprintf('arrival %.2f m/s', r.flyVms); end
    if r.massErr > 1e-6,      bad{end+1} = sprintf('mass law %.1e', r.massErr); end

    % WITNESS
    [okW, za] = cap(pool, 300, @pumpkyn.cr3bp.tfMin, 1, rv0(1:6)', rvf(1:6)', z8, Tnd, cnd, mu);
    if okW && isnumeric(za) && numel(za) == 8 && all(isfinite(za))
        r.dz = norm(za(:) - z8(:));
        [okWF, ~, Ya] = cap(pool, 300, @pumpkyn.cr3bp.tfMinProp, 2, za(8), [rv0(1:6); 1; za(1:7)], Tnd, cnd, mu);
        if okWF, r.witKm = norm(Ya(end,1:3) - rvf(1:3)')*lStar; end
    end
    if ~(isfinite(r.dz) && r.dz <= tolDz), bad{end+1} = sprintf('witness |dz| %.1e', r.dz); end
    if ~(isfinite(r.witKm) && r.witKm < gateKm), bad{end+1} = sprintf('witness flight %.1f km', r.witKm); end

    % SECOND ORDER, recomputed
    seed = seed_from_z8(z8, rv0(1:6), 24, Tnd, cnd, mu);
    % FAIL CLOSED. Every check below used to read `isfinite(x) && x ~= ...`,
    % so a re-polish that timed out (x = NaN) was a CLEAN row: the audit
    % passed exactly the entries it had failed to examine (FINDINGS 78). A
    % result that is missing, malformed or not a pass is now a BAD row.
    [okP, zp, it] = cap(pool, 900, @ms_tfmin, 2, rv0(1:6), rvf(1:6), seed, Tnd, cnd, mu, ...
                           struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));
    if ~(okP && isnumeric(zp) && numel(zp) == 8 && all(isfinite(zp)))
        bad{end+1} = 're-polish failed or timed out: no second-order verdict';
    else
        % the verdict below belongs to the POLISHED root; it speaks for the
        % stored one only if the polish did not leave it
        r.dzPolish = norm(zp(:) - z8(:));
        if ~(r.dzPolish <= tolMove), bad{end+1} = sprintf('polished root moved %.1e from the stored one (limit %g)', r.dzPolish, tolMove); end
        if isstruct(it) && isfield(it, 'conj') && isfield(it.conj, 'pass')
            [okc, cv] = scalar_verdict(it.conj.pass);
            if okc, r.conjNow = cv; end
        end
        if ~isfinite(r.conjNow),  bad{end+1} = 're-polish returned no conjugate verdict';
        elseif r.conjNow ~= 1,    bad{end+1} = sprintf('conjugate verdict %g', r.conjNow);
        end
        if isfinite(r.conjStored) && isfinite(r.conjNow) && r.conjStored ~= r.conjNow
            bad{end+1} = sprintf('conj stored %g, recomputed %g', r.conjStored, r.conjNow);
        end
    end
    [okG, g] = cap(pool, 900, @mintime_hypothesis_gates, 1, z8, rv0(1:6), Tnd, cnd, mu, struct());
    if ~(okG && isstruct(g) && isscalar(g) && all(isfield(g, {'dimS', 'minLamVBound', 'minQmtBound', 'h6Ok'})))
        bad{end+1} = 'hypothesis gates failed or timed out';
    else
        % VALIDATE, THEN COMPARE (as the certifier does). `x > floor` is true for
        % Inf, logical(2) is true, and logical(NaN) THROWS -- so each value must
        % be a real finite scalar before it is read as a measurement.
        if ~realScalar(g.minLamVBound) || ~(g.minLamVBound > hypFloor), bad{end+1} = sprintf('H2: lower bound of min|lam_v| over the arc is %s, not a finite value above %g', showValue(g.minLamVBound), hypFloor); end
        if ~realScalar(g.minQmtBound)  || ~(g.minQmtBound  > hypFloor), bad{end+1} = sprintf('H3: lower bound of min Q_mt over the arc is %s, not a finite value above %g', showValue(g.minQmtBound), hypFloor); end
        if ~realScalar(g.dimS) || g.dimS ~= 1
            bad{end+1} = sprintf('dim S = %s, not 1', showValue(g.dimS));
        else
            r.dimSNow = g.dimS;
        end
        if ~(realScalar(g.h6Ok) && g.h6Ok == 1), bad{end+1} = sprintf('H6 does not hold (h6Ok = %s)', showValue(g.h6Ok)); end
        if isfinite(r.dimSStored) && isfinite(r.dimSNow) && r.dimSStored ~= r.dimSNow
            bad{end+1} = sprintf('dim S stored %g, recomputed %g', r.dimSStored, r.dimSNow);
        end
    end

    r.ok = isempty(bad);
    if ~r.ok, r.problem = strjoin(bad, '; ');  problems{end+1} = sprintf('(%d,%d): %s', r.iD, r.iA, r.problem); end %#ok<AGROW>
    rows(end+1) = r; %#ok<AGROW>
    fprintf('AUDIT (%2d,%2d) sD %.4f sA %.4f  t_f %.6f  fly %.3f km / %.3f m/s  dz %.1e  wit %.3f km  conj %g/%g  %s\n', ...
        r.iD, r.iA, r.sD, r.sA, r.tfZ8, r.flyKm, r.flyVms, r.dz, r.witKm, r.conjStored, r.conjNow, ...
        ternS(r.ok, 'OK', ['BAD: ' r.problem]));
end

A = struct('rows', rows, 'nOk', nnz([rows.ok]), 'nBad', nnz(~[rows.ok]), ...
           'problems', {problems}, 'catMat', catMat, 'when', datestr(now));
fprintf('AUDIT SUMMARY: %d entries audited, %d OK, %d BAD\n', numel(rows), A.nOk, A.nBad);
if ~isempty(d('out', '')), save(d('out', ''), 'A'); end
end

function varargout = capped(pool, capSec, fh, nout, varargin)
% CAPPED  One external call under a hard cap (direct when no pool).
% INPUTS: pool; capSec; fh; nout; varargin.  OUTPUTS: ok, then nout outputs.
varargout = cell(1, nout + 1);
if isempty(pool)
    try
        [varargout{2:nout+1}] = feval(fh, varargin{:});  varargout{1} = true;
    catch
        varargout{1} = false;
    end
    return
end
[varargout{1}, varargout{2:nout+1}] = run_capped(pool, fh, nout, capSec, varargin{:});
end

function tf = realScalar(x)
% REALSCALAR  Is x one real, finite number (or one logical)? The test a gate
% value must pass BEFORE it is compared with anything.  INPUTS: x.
% OUTPUTS: tf.
tf = (isnumeric(x) || islogical(x)) && isscalar(x) && isreal(x) && isfinite(double(x));
end

function t = showValue(x)
% SHOWVALUE  A gate value as text for a BAD row, whatever its shape.
% INPUTS: x.  OUTPUTS: t char.
if (isnumeric(x) || islogical(x)) && numel(x) <= 4, t = mat2str(double(x), 4); else, t = sprintf('<%s %s>', class(x), mat2str(size(x))); end
end

function s = ternS(c, a, b)
% TERNS  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
