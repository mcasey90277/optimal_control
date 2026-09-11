function R = report_optimality(T, opts)
%% Purpose:
%
%   The optimality report for ONE transfer, in the two sections the theory
%   actually has.
%
%   SECTION 1 -- NECESSARY (Pontryagin, first order). The multiple-shooting
%   residual (the costate equations, terminal matching and transversality),
%   the Hamiltonian along the arc, the flown arrival, and an independent
%   re-solve. A trajectory can satisfy every line here and still not be a
%   minimizer: on this very problem, 12 of 14 candidates did exactly that.
%
%   SECTION 2 -- SUFFICIENCY HYPOTHESES (Bonnard-Caillau-Trelat). Normality
%   (no abnormal lift of the same trajectory), strengthened Legendre, the
%   all-burn switching function, and no conjugate time. With SECTION 1 these
%   give a strict strong local minimizer among trajectories with the same
%   endpoints -- numerically, at the sampled times.
%
%   Section 1 is NOT redundant given section 2: being an extremal is a
%   HYPOTHESIS of the theorem, not a consequence, and the conjugate test is
%   computed along the extremal -- off one, its determinant means nothing.
%   So a failed section 1 suppresses the claim outright.
%
%% Inputs:
%
%  T                        struct                  a certify_root output
%  opts                     struct (optional)
%   .quiet [false] .tolR [3e-11] .gateKm [100] .gateVms [10] .tolDz [1e-6]
%   .tolH [1e-6] Hamiltonian residual along the arc
%
%% Outputs:
%
%  R                        struct                  .necessary .sufficient
%                                                   (each requires its lines
%                                                   to be CHECKED and to
%                                                   pass -- "all checked
%                                                   passed" is vacuously
%                                                   true when nothing was)
%                                                   .necessaryComplete
%                                                   .sufficientComplete
%                                                   .claim (both) .verdict
%                                                   (text) .lines (cellstr)
%                                                   .rows (name/value/ok)
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
tolR = d('tolR', 3e-11);  gateKm = d('gateKm', 100);  gateVms = d('gateVms', 10);
tolDz = d('tolDz', 1e-6);  tolH = d('tolH', 1e-6);
g = [];  if isfield(T, 'g') && isstruct(T.g) && ~isempty(T.g), g = T.g; end
lines = {};  rows = struct('name', {}, 'value', {}, 'ok', {}, 'section', {});

add = @(nm, val, okv, sec) deal(nm, val, okv, sec);   %#ok<NASGU>
    function push(nm, val, okv, sec, fmt)
        rows(end+1) = struct('name', nm, 'value', val, 'ok', okv, 'section', sec); %#ok<AGROW>
        if isnan(okv)
            st = 'NOT CHECKED';
        elseif okv
            st = 'PASS';
        else
            st = 'FAIL';
        end
        if isnan(val), vs = '     --   '; else, vs = sprintf(fmt, val); end
        lines{end+1} = sprintf('  %-46s %s   %s', nm, vs, st); %#ok<AGROW>
    end

% ---- section 1: necessary ------------------------------------------------
lines{end+1} = 'NECESSARY (Pontryagin, first order)';
push('BVP residual (costates, terminal, transversality)', gv(T, 'normR'), ...
     tri(gv(T,'normR'), @(x) x <= max(tolR, 1e-8)), 1, '%9.2e');
push('Hamiltonian along the arc  |lambda.f + 1|', gg(g, 'Hresid'), ...
     tri(gg(g,'Hresid'), @(x) x <= tolH), 1, '%9.2e');
push('flown arrival, position [km]', gv(T, 'flyKm'), ...
     tri(gv(T,'flyKm'), @(x) x < gateKm), 1, '%9.4f');
push('flown arrival, velocity [m/s]', gv(T, 'flyVms'), ...
     tri(gv(T,'flyVms'), @(x) x < gateVms), 1, '%9.4f');
push('independent solver agreement  |dz|', gv(T, 'dz'), ...
     tri(gv(T,'dz'), @(x) x <= tolDz), 1, '%9.2e');
push('that solver''s own flight [km]', gv(T, 'flyKmWitness'), ...
     tri(gv(T,'flyKmWitness'), @(x) x < gateKm), 1, '%9.4f');

% ---- section 2: sufficiency hypotheses -----------------------------------
lines{end+1} = '';
lines{end+1} = 'SUFFICIENCY HYPOTHESES (Bonnard-Caillau-Trelat)';
push('H1 normality: dim S (1 = no abnormal lift)', gg(g, 'dimS'), ...
     tri(gg(g,'dimS'), @(x) x == 1), 2, '%9.0f');
push('H2 strengthened Legendre: min |lambda_v| > 0', gg(g, 'minLamV'), ...
     tri(gg(g,'minLamV'), @(x) x > 0), 2, '%9.3e');
push('H3 all-burn is the PMP control: min Q > 0', gg(g, 'minQmt'), ...
     tri(gg(g,'minQmt'), @(x) x > 0), 2, '%9.3e');
push('H5 no conjugate time in (0, t_f]', gv(T, 'conj'), ...
     tri(gv(T,'conj'), @(x) x == 1), 2, '%9.0f');

% A section holds only if every one of its lines was CHECKED and passed.
% Vacuous truth is the failure mode here: "all checked lines passed" is
% trivially true when nothing was checked, and that is exactly the sentence
% a report must never be able to produce.
s1 = [rows.section] == 1;  s2 = [rows.section] == 2;
o1 = okOf(rows(s1));  o2 = okOf(rows(s2));
R = struct();
R.necessaryComplete  = ~any(isnan(o1));
R.sufficientComplete = ~any(isnan(o2));
R.necessary  = R.necessaryComplete  && all(o1 == 1);
R.sufficient = R.sufficientComplete && all(o2 == 1);
R.claim      = R.necessary && R.sufficient;
if R.claim
    R.verdict = ['strict strong local minimizer among trajectories with the same ' ...
                 'endpoints -- numerically certified at the sampled times'];
elseif ~R.necessaryComplete
    R.verdict = ['NOT CHECKED: a necessary condition was never evaluated, so nothing ' ...
                 'is claimed (an unchecked line is not a passed line)'];
elseif ~all(o1 == 1)
    R.verdict = ['NOT an extremal to tolerance: the second-order test is meaningless ' ...
                 'off an extremal, so no minimality is claimed'];
elseif ~R.sufficientComplete
    R.verdict = 'an extremal; sufficiency NOT CHECKED (missing diagnostics)';
else
    R.verdict = 'an extremal, but a sufficiency hypothesis fails: no minimality claimed';
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
