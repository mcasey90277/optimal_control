function R = report_optimality(T, opts)
%% Purpose:
%
%   The optimality report for ONE transfer, in the three groups the theory
%   and the instruments actually have.
%
%   NECESSARY (Pontryagin, first order): the multiple-shooting residual (the
%   costate equations, terminal matching and transversality), the
%   Hamiltonian along the arc, the flown arrival, and -- when the caller
%   supplies them -- the adjoint-equation check (N5) and the exact
%   minimum-principle gap of the applied control (N6). A trajectory can
%   satisfy every line here and still not be a minimizer: on this problem,
%   12 of 14 candidates did exactly that.
%
%   SUFFICIENCY HYPOTHESES (Bonnard-Caillau-Trelat): normality (no
%   abnormal lift of the same trajectory), strengthened Legendre, the
%   all-burn switching function, and no conjugate time. With the necessary
%   lines these give a strict strong local minimizer among trajectories
%   with the same endpoints. H6 (lambda_m(0) < c/T with margin) is listed
%   here as the VALIDITY condition of the reduced conjugate instrument: a
%   conjugate verdict without it is not interpretable.
%
%   CROSS-CHECKS: the independent solver's agreement and its own flight.
%   These are guards against OUR implementation, not conditions of the
%   theory. They are kept out of the PMP conjunction -- a failed cross-check
%   never reads as "not an extremal" -- but they still block the claim.
%
%   The necessary lines are NOT redundant given sufficiency: being an
%   extremal is a HYPOTHESIS of the theorem, not a consequence, and the
%   conjugate test is computed along the extremal -- off one, its
%   determinant means nothing.
%
%   Every line has FOUR states: PASS, FAIL, NOT CHECKED (the quantity was
%   never supplied) and UNRESOLVED (the instrument ran and could not decide,
%   e.g. a conjugate ENDPOINT verdict). Only PASS counts toward the claim;
%   the verdict text says which of the other three blocked it.
%
%% Inputs:
%
%  T                        struct                  a certify_root output;
%                                                   optional .adjErr .dirGap
%                                                   .lamMf .h6Margin
%                                                   .conjVerdict
%                                                   .flyVmsWitness
%  opts                     struct (optional)
%   .quiet [false] .tolR [3e-11] .gateKm [100] .gateVms [10] .tolDz [1e-6]
%   .tolH [1e-6] Hamiltonian residual, .tolAdj [1e-7] adjoint relative
%   error, .tolGap [1e-12] minimum-principle gap, .tolLamMf [1e-6]
%   transversality, .h6MarginMin [1]
%
%% Outputs:
%
%  R                        struct                  .necessary .sufficient
%                                                   .crossChecks (each
%                                                   requires ALL its lines
%                                                   to be PASS -- "all
%                                                   checked lines passed"
%                                                   is vacuously true when
%                                                   nothing was)
%                                                   .necessaryComplete
%                                                   .sufficientComplete
%                                                   .unresolved .claim
%                                                   .verdict (text) .lines
%                                                   (cellstr) .rows
%                                                   (id/name/value/ok/
%                                                   section)
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
tolR = d('tolR', 3e-11);  gateKm = d('gateKm', 100);  gateVms = d('gateVms', 10);
tolDz = d('tolDz', 1e-6);  tolH = d('tolH', 1e-6);  tolAdj = d('tolAdj', 1e-7);
tolGap = d('tolGap', 1e-12);  tolLamMf = d('tolLamMf', 1e-6);
h6MarginMin = d('h6MarginMin', 1);
g = [];  if isfield(T, 'g') && isstruct(T.g) && ~isempty(T.g), g = T.g; end
lines = {};
rows = struct('id', {}, 'name', {}, 'value', {}, 'ok', {}, 'section', {});

% ok codes: 1 PASS, 0 FAIL, NaN NOT CHECKED, -1 UNRESOLVED
    function push(id, nm, val, okv, sec, fmt)
        rows(end+1) = struct('id', id, 'name', nm, 'value', val, 'ok', okv, 'section', sec); %#ok<AGROW>
        if isnan(okv),      st = 'NOT CHECKED';
        elseif okv == -1,   st = 'UNRESOLVED';
        elseif okv == 1,    st = 'PASS';
        else,               st = 'FAIL';
        end
        if isnan(val), vs = '     --   '; else, vs = sprintf(fmt, val); end
        lines{end+1} = sprintf('  %-3s %-44s %s   %s', id, nm, vs, st); %#ok<AGROW>
    end

% ---- NECESSARY --------------------------------------------------------
lines{end+1} = 'NECESSARY (Pontryagin, first order)';
push('N1', 'BVP residual (costates, terminal, transversality)', gv(T, 'normR'), ...
     tri(gv(T,'normR'), @(x) x <= max(tolR, 1e-8)), 1, '%9.2e');
push('N2', 'Hamiltonian along the arc  |lambda.f + 1|', firstOf(gv(T, 'Hresid'), gg(g, 'Hresid')), ...
     tri(firstOf(gv(T, 'Hresid'), gg(g, 'Hresid')), @(x) x <= tolH), 1, '%9.2e');
push('N3', 'flown arrival, position [km]', gv(T, 'flyKm'), ...
     tri(gv(T,'flyKm'), @(x) x < gateKm), 1, '%9.4f');
push('N3', 'flown arrival, velocity [m/s]', gv(T, 'flyVms'), ...
     tri(gv(T,'flyVms'), @(x) x < gateVms), 1, '%9.4f');
push('N4', 'transversality |lambda_m(t_f)|', gv(T, 'lamMf'), ...
     tri(gv(T,'lamMf'), @(x) x <= tolLamMf), 1, '%9.2e');
push('N5', 'adjoint equations, relative error', gv(T, 'adjErr'), ...
     tri(gv(T,'adjErr'), @(x) x <= tolAdj), 1, '%9.2e');
push('N6', 'minimum principle, gap of the applied control', gv(T, 'dirGap'), ...
     tri(gv(T,'dirGap'), @(x) x <= tolGap), 1, '%9.2e');

% ---- SUFFICIENCY HYPOTHESES ------------------------------------------
lines{end+1} = '';
lines{end+1} = 'SUFFICIENCY HYPOTHESES (Bonnard-Caillau-Trelat)';
push('S1', 'strengthened Legendre: min |lambda_v| > 0', gg(g, 'minLamV'), ...
     tri(gg(g,'minLamV'), @(x) x > 0), 2, '%9.3e');
push('S2', 'strict bang: min Q_mt > 0', gg(g, 'minQmt'), ...
     tri(gg(g,'minQmt'), @(x) x > 0), 2, '%9.3e');
push('S3', 'normality: dim S (1 = no abnormal lift)', gg(g, 'dimS'), ...
     tri(gg(g,'dimS'), @(x) x == 1), 2, '%9.0f');
conjOk = tri(gv(T,'conj'), @(x) x == 1);
if conjOk == 0 && isfield(T, 'conjVerdict') && strcmpi(T.conjVerdict, 'ENDPOINT')
    conjOk = -1;                       % the instrument ran and could not decide
end
push('S4', 'no conjugate time in (0, t_f]', gv(T, 'conj'), conjOk, 2, '%9.0f');
h6 = firstOf(gv(T, 'h6Margin'), gg(g, 'h6Margin'));
h6ok = tri(h6, @(x) x > h6MarginMin);            % strict, as h6_margin itself
if h6ok == 1 && isfield(g, 'h6Ok') && ~isempty(g.h6Ok) && ~g.h6Ok
    h6ok = 0;                                    % the helper's clearance verdict governs
end
push('S5', 'H6 instrument validity: (c/T)/lambda_m(0) margin', h6, h6ok, 2, '%9.2f');

% ---- CROSS-CHECKS -----------------------------------------------------
lines{end+1} = '';
lines{end+1} = 'CROSS-CHECKS (implementation guards, not PMP conditions)';
push('X1', 'independent solver agreement  |dz|', gv(T, 'dz'), ...
     tri(gv(T,'dz'), @(x) x <= tolDz), 3, '%9.2e');
push('X2', 'that solver''s own flight, position [km]', gv(T, 'flyKmWitness'), ...
     tri(gv(T,'flyKmWitness'), @(x) x < gateKm), 3, '%9.4f');
push('X2', 'that solver''s own flight, velocity [m/s]', gv(T, 'flyVmsWitness'), ...
     tri(gv(T,'flyVmsWitness'), @(x) x < gateVms), 3, '%9.4f');

% A group holds only if every one of its lines was CHECKED and PASSED.
% Vacuous truth is the failure mode here: "all checked lines passed" is
% trivially true when nothing was checked, and that is exactly the sentence
% a report must never be able to produce.
o1 = okOf(rows([rows.section] == 1));
o2 = okOf(rows([rows.section] == 2));
o3 = okOf(rows([rows.section] == 3));
R = struct();
R.necessaryComplete  = ~any(isnan(o1));
R.sufficientComplete = ~any(isnan(o2));
R.unresolved         = any(o2 == -1);
R.necessary   = R.necessaryComplete  && all(o1 == 1);
R.sufficient  = R.sufficientComplete && all(o2 == 1);
R.crossChecks = ~any(isnan(o3)) && all(o3 == 1);
R.claim       = R.necessary && R.sufficient && R.crossChecks;
if R.claim
    R.verdict = ['All required numerical checks passed. If the stated hypotheses hold ' ...
                 'exactly, this arc is a strict strong local minimizer among trajectories ' ...
                 'with the same endpoints. The evidence is numerical and sampled.'];
elseif ~R.necessaryComplete
    R.verdict = ['NOT CHECKED: a necessary condition was never evaluated, so nothing ' ...
                 'is claimed (an unchecked line is not a passed line)'];
elseif ~all(o1 == 1)
    R.verdict = ['NOT an extremal to tolerance: the second-order test is meaningless ' ...
                 'off an extremal, so no minimality is claimed'];
elseif R.unresolved
    R.verdict = ['an extremal; the conjugate test is UNRESOLVED (ENDPOINT verdict): ' ...
                 'no minimality is claimed and none is refuted'];
elseif ~R.sufficientComplete
    R.verdict = 'an extremal; sufficiency NOT CHECKED (a hypothesis or H6 was never evaluated)';
elseif ~all(o2 == 1)
    R.verdict = 'an extremal, but a sufficiency hypothesis fails: no minimality claimed';
else
    R.verdict = ['every PMP and sufficiency line passed, but a CROSS-CHECK failed: the ' ...
                 'implementation is in doubt, so no claim is made'];
end
lines{end+1} = '';
lines{end+1} = sprintf('VERDICT: %s', R.verdict);
R.lines = lines(:);  R.rows = rows;

if ~d('quiet', false)
    fprintf('\n=========== OPTIMALITY REPORT ===========\n');
    if isfield(T, 'sD') && isfield(T, 'sA')
        fprintf('  departure phase %.4f   arrival phase %.4f\n', T.sD, T.sA);
    end
    if isfield(T, 'tfDays')
        fprintf('  t_f = %.4f d   dV = %.4f km/s   propellant %.2f kg\n\n', ...
                T.tfDays, gv(T,'dvKms'), gv(T,'propellantKg'));
    end
    fprintf('%s\n', lines{:});
    fprintf('=========================================\n\n');
end
end

function v = gv(T, f)
% GV  Numeric field or NaN.  INPUTS: T; f.  OUTPUTS: v.
v = NaN;
if isfield(T, f)
    [okv, x] = scalar_verdict(T.(f));
    if okv, v = x; end
end
end

function v = gg(g, f)
% GG  Numeric gate field or NaN.  INPUTS: g; f.  OUTPUTS: v.
v = NaN;
if ~isempty(g) && isstruct(g) && isfield(g, f)
    [okv, x] = scalar_verdict(g.(f));
    if okv, v = x; end
end
end

function v = firstOf(a, b)
% FIRSTOF  The first non-NaN of two candidates.  INPUTS: a; b.  OUTPUTS: v.
if ~isnan(a), v = a; else, v = b; end
end

function o = tri(v, test)
% TRI  PASS / FAIL / NOT CHECKED as 1 / 0 / NaN. A missing value is NEVER a
% pass.  INPUTS: v; test.  OUTPUTS: o.
if isnan(v), o = NaN; else, o = double(test(v)); end
end

function o = okOf(rr)
% OKOF  The ok column.  INPUTS: rr.  OUTPUTS: o.
if isempty(rr), o = []; else, o = [rr.ok]; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
