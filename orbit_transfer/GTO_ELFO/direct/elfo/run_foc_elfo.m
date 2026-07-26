function [rep, info] = run_foc_elfo(kind, matFile)
% RUN_FOC_ELFO  Standard first-order optimality report for a certified ELFO
% row -- the FIRST first-order optimality gate for the GTO->ELFO campaign.
%
% Thin wiring driver (Task 9, ELFO precinct of the FOC-gate campaign): loads a
% certified ELFO row, warm re-solves it AT the saved primal with the CasADi
% model attached (casadi_energy_freetf.m / casadi_mintime_freetf.m, Task 9's
% opts.returnModel hook), guards on certified quantities, runs the generic
% AD-based FOC/KKT gate (foc_check.m) with the IPOPT-inertia 2nd-order verdict
% attached (foc_ipopt_inertia.m), and prints the standard report (foc_report.m).
%
% Two kinds, selecting solver + manifest:
%   'fuel'    - casadi_energy_freetf.m / foc_manifest('elfo_fuel')  (nu=4,
%               throttle row present -- sign-law/singular/Sdot checks run).
%               Handles BOTH a pure energy seed (epsilon=1, e.g.
%               energy_elfo_freetf.mat, no top-level .epsilon field -> default
%               1) and a sharpened/bang-bang fuel row (epsilon from the file,
%               e.g. minfuel_ELFO_tf*_minEps0.mat). The bang-bang sign law is
%               only EXACT at epsilon=0 (Bertrand-Epenoy smoothing at
%               epsilon>0 makes S<0<=>full-thrust an approximate law); at
%               epsilon>0 rep.signPct is printed with an explicit caveat and
%               excluded from the reader's mental "pass" bar (foc_check's own
%               rep.pass still folds it in -- see rep.epsCaveat / the printed
%               NOTE for the honest read).
%   'mintime' - casadi_mintime_freetf.m / foc_manifest('elfo_mintime')  (nu=3,
%               thrRow=[] -- hard all-burn, sign-law/singular/Sdot checks
%               correctly SKIP). The free-t_f dual-form horizon condition
%               (G4) is what rep.lamTimeCoV / rep.lamTimeEnd report here.
%
% INPUTS:
%   kind    - 'fuel' | 'mintime' [char]
%   matFile - certified row (.mat) [char, default: results/energy_elfo_freetf.mat
%             for 'fuel', results/mintime_elfo.mat for 'mintime']. 'fuel' file
%             must carry (at least): X[9x(N+1)], U[4x(N+1)], sigma, rv0, rvf,
%             tauf0, moonZone, pSund, qSund (and optionally epsilon -- absent
%             means an epsilon=1 energy seed). 'mintime' file must carry:
%             X[9x(N+1)], U[3 or 4 x(N+1)], sigma, rv0, rvf, tauf0, tf,
%             moonZone, pSund, qSund (elfo_run_one / gen_elfo_minfuel /
%             gen_elfo_mintime save-field lists).
%
% OUTPUTS:
%   rep  - foc_check report struct, with .ipopt (2nd-order verdict) attached,
%          and (fuel, epsilon>0 only) .epsCaveat [char] [struct]
%   info - warm re-solve bookkeeping: .drift .ipoptStatus .maxDefect .epsilon
%          (fuel: .mfSaved .mfResolved) (mintime: .tfSaved .tfResolved)
%          .improved [struct]
%
% REFERENCES:
%   [1] verify_common/foc_check.m, foc_report.m, foc_ipopt_inertia.m.
%   [2] casadi_energy_freetf.m / casadi_mintime_freetf.m (Task 9 returnModel/
%       creg + regHistory hook).
%   [3] GTO_tulip/direct/sundman_minfuel/run_foc_tulip.m (Task 8 wrapper shape
%       this mirrors: warm re-solve at saved primal + certified-quantity
%       guard + foc_check + rep.ipopt + foc_report). The tulip driver's
%       LS-reconstructed-costate cross-check (certify_minfuel_pmp.m) is
%       SINGLE-primary fixed-t_f machinery that does not apply to this
%       two-primary free-t_f model (run_elfo_minfuel.m stage 3 note) and is
%       not ported here -- foc_check's raw-dual costate is the only source.
%   [4] gen_elfo_minfuel.m / gen_elfo_mintime.m (the artifacts' own solve
%       configs this warm re-solve reuses: moonZone/pSund/qSund/cBox/
%       tfCapMult/tfTarget).

here = fileparts(mfilename('fullpath'));
addpath(here);  setup_paths();               % sundman engine + cr3bp_common
vcDir = fullfile(here, '..', '..', '..', 'verify_common');
addpath(vcDir);
setup_verify_common();

assert(ischar(kind) || isstring(kind), 'run_foc_elfo: kind must be char/string');
kind = lower(char(kind));
assert(any(strcmp(kind, {'fuel','mintime'})), ...
    'run_foc_elfo: kind must be ''fuel'' or ''mintime'' (got ''%s'')', kind);

if nargin < 2 || isempty(matFile)
    if strcmp(kind, 'fuel')
        matFile = fullfile(here, 'results', 'energy_elfo_freetf.mat');
    else
        matFile = fullfile(here, 'results', 'mintime_elfo.mat');
    end
end

p = cr3bp_lt_params(0.025, 15, 2100);   % nominal 25 mN rung (matches the certified campaign)
S = load(matFile);
feasTol = 1e-8;

switch kind
% =============================================================================
case 'fuel'
    man = foc_manifest('elfo_fuel');
    epsilon = 1;  if isfield(S,'epsilon') && ~isempty(S.epsilon), epsilon = S.epsilon; end
    tf0 = S.X(8,end);

    % Warm re-solve AT the saved primal, PINNED to the same t_f (tfTarget --
    % well-posed, matches gen_elfo_minfuel's own step_solve convention) and
    % the artifact's own two-primary clock / cBox / tfCapMult config.
    sopts = struct('returnModel', true, 'moonZone', S.moonZone, 'pSund', S.pSund, ...
        'qSund', S.qSund, 'epsilon', epsilon, 'tfTarget', tf0, 'cBox', [0.15 6], ...
        'tfCapMult', 6, 'warmTight', true, 'maxIter', 2000);
    out = casadi_energy_freetf(S.sigma, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
        S.X, S.U, S.tauf0, sopts);

    % --- certified-quantity guard (shared: verify_common/certified_guard.m) --
    % Success class + machine-tight defect + ONE-SIDED final-mass check.
    % Extracted 2026-07-26 from four identical copies across four campaigns.
    mfSaved = S.X(7,end);
    info = certified_guard( ...
        struct('success', out.success, 'ipoptStatus', out.ipoptStatus, ...
               'maxDefect', out.maxDefect, 'value', out.mf), ...
        struct('caller','run_foc_elfo', 'label', matFile, 'saved', mfSaved, ...
               'name','m_f', 'errName','mass', 'better','higher', ...
               'feasTol', feasTol, 'tol', 1e-6));
    info.drift      = max(abs(out.X(:) - S.X(:)));
    info.epsilon    = epsilon;
    info.mfSaved    = mfSaved;
    info.mfResolved = out.mf;

    % Pass the epsilon this row was re-solved at: foc_check skips the
    % bang-bang family (sign law / singular arc / regular switching) when
    % the throttle cost is not affine, instead of folding meaningless
    % numbers into the advisory verdict (external review, 2026-07-25).
    rep = foc_check(out, S.sigma, man, struct('eps', epsilon));
    rep.ipopt = foc_ipopt_inertia(getfield_default(out, 'regHistory', []));

    if epsilon > 0
        rep.epsCaveat = sprintf(['epsilon=%.3g leg: the throttle sign law (S<0 <=> full ' ...
            'thrust) is only EXACT at epsilon=0 (bang-bang). At epsilon=%.3g the Bertrand-' ...
            'Epenoy smoothing makes the throttle a continuous function of S, so ' ...
            'rep.signPct (%.2f%%) is REPORT-ONLY here, not a pass/fail discriminator.'], ...
            epsilon, epsilon, rep.signPct);
        fprintf('\n  [run_foc_elfo] NOTE: %s\n', rep.epsCaveat);
    end

% =============================================================================
case 'mintime'
    man = foc_manifest('elfo_mintime');
    tfSaved = S.X(8,end);  if isfield(S,'tf') && ~isempty(S.tf), tfSaved = S.tf; end

    sopts = struct('returnModel', true, 'moonZone', S.moonZone, 'pSund', S.pSund, ...
        'qSund', S.qSund, 'cBox', [0.10 8], 'tfCapMult', 4, 'warmTight', true, 'maxIter', 3000);
    out = casadi_mintime_freetf(S.sigma, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
        S.X, S.U, S.tauf0, sopts);

    % --- certified-quantity guard (shared: verify_common/certified_guard.m) --
    % Min-time has no "one-sided mass" analog -- it does not optimize mass at
    % all. The certified quantity here is t_f, and the one-sided test INVERTS:
    % an EQUAL or LOWER t_f is fine (that is what min-time is for); a HIGHER
    % one means the recovered point is a worse anchor and must be refused.
    % That inversion is what spec.better='lower' selects, and it is pinned by
    % verify_common/tests/test_certified_guard.m.
    info = certified_guard( ...
        struct('success', out.success, 'ipoptStatus', out.ipoptStatus, ...
               'maxDefect', out.maxDefect, 'value', out.tf), ...
        struct('caller','run_foc_elfo', 'label', matFile, 'saved', tfSaved, ...
               'name','t_f', 'errName','tf', 'better','lower', ...
               'feasTol', feasTol, 'tol', 1e-6, 'units', ' ND'));
    info.drift      = max(abs(out.X(:) - S.X(:)));
    info.tfSaved    = tfSaved;
    info.tfResolved = out.tf;

    rep = foc_check(out, S.sigma, man, struct());
    rep.ipopt = foc_ipopt_inertia(getfield_default(out, 'regHistory', []));
end

[~, tag] = fileparts(matFile);
resDir = fullfile(here, 'results');
foc_report(rep, tag, resDir);
end

% =============================================================================
function v = getfield_default(s, f, dflt)
% GETFIELD_DEFAULT  s.(f) if present and nonempty, else dflt.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
