function C = certify_dro_mintime(o, p, Tmax, c, opts)
% CERTIFY_DRO_MINTIME  Pass/fail gate for a direct DRO->tulip min-time solution.
%
% WHY THIS EXISTS. Every solve in this campaign reported a defect of 1e-14 and
% every one of them was wrong. The defect is a statement about the
% discretization; it is not evidence that the trajectory is physical. This gate
% is the standard a solve must meet before its t_f is quoted as an answer.
%
% THE TWO GATES THAT MATTER, and both are independent of the solver:
%
%   G1  LOCAL ACCURACY. The worst continuous-time local error must be below a
%       position tolerance stated in KILOMETRES, not in solver units. Default
%       1 km. Measured by re-integrating each interval at 1e-11/1e-13 -- see
%       dro_residual.
%
%   G1b GLOBAL ACCURACY, and this is the one a mission cares about. G1 restarts
%       at every node, so it is a max over LOCAL errors and could in principle
%       understate what happens when the errors accumulate. G1b instead flies
%       the reconstructed control ONCE, from the initial node to the end, and
%       asks where the spacecraft actually arrives. Default 1 km.
%
%       Measured on the certified N = 1600 Hermite-Simpson solution: max local
%       0.322 km, SUM of local 1.746 km, global single-shot 0.700 km -- i.e.
%       the errors partially CANCEL rather than accumulate, and the terminal
%       position error is 0.0 km with 0.002 m/s of velocity error. Worth
%       stating because the opposite was the worry: over 1600 intervals in an
%       unstable three-body field, accumulation was the expected failure mode
%       and it did not occur.
%
%   G2  AGREEMENT. This is the one problem in the catalog where an independently
%       converged INDIRECT solution exists (pumpkyn.cr3bp.tfMin, t_f = 4.015243).
%       A working direct solver must reproduce it. Default tolerance 1e-4
%       relative, i.e. agreement to four digits.
%
% G2 is the sharper test and it is the reason DRO->tulip was chosen to start
% with. A direct method that cannot reproduce a known answer on the one problem
% where the answer is known should not be trusted on the problems where it is
% not.
%
% THE SUPPORTING CHECKS (G3-G6) are the usual transcription hygiene: defect,
% control unit norm, terminal boundary error, throttle saturation, the lifted-
% time spread, and -- when a floor was imposed -- whether the trajectory honours
% it BETWEEN nodes and not merely at them. They are necessary, not sufficient,
% and they are reported separately so that a solve which passes all of them and
% still fails G1 is visible for what it is.
%
% INPUTS:
%   o     - solution struct from casadi_mintime_dro
%   p     - params from dro_tulip_endpoints (.muStar .lStar .tStar)
%   Tmax  - ND thrust acceleration at unit mass fraction [scalar]
%   c     - ND exhaust velocity [scalar]
%   opts  - struct (optional):
%           .posTolKm  [1.0]      G1 tolerance, km
%           .globTolKm [1.0]      G1b tolerance, km. [] skips G1b (it costs one
%                                 full-trajectory propagation).
%           .tfRef     [4.0152430] G2 reference t_f, ND. [] skips G2.
%           .tfRelTol  [1e-4]     G2 relative tolerance
%           .defectTol [1e-9]     G3
%           .rMoonKm   [1737.4]
%           .verbose   [true]
%
% OUTPUTS:
%   C - struct: .pass (G1 && G2), .passAll (every gate), .gates (struct array
%       with .id .name .value .tol .pass .units), .resid (from dro_residual),
%       .worstAltKm, .tfErrRel
%
% REFERENCES:
%   [1] orbit_transfer/DRO_tulip/FINDINGS.md
%   [2] orbit_transfer/OPTIMALITY_CERTIFICATION.md -- the repo-wide register
%       this gate reports into.

if nargin < 5, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
posTolKm  = d('posTolKm', 1.0);
globTolKm = d('globTolKm', 1.0);
tfRef     = d('tfRef', 4.0152430);
tfRelTol  = d('tfRelTol', 1e-4);
defectTol = d('defectTol', 1e-9);
rMoonKm   = d('rMoonKm', 1737.4);
verbose   = d('verbose', true);

lStar = p.lStar;  mu = p.muStar;

% G2 compares against pumpkyn.cr3bp.tfMin, which has NO path-constraint
% capability. Once a floor is imposed the constrained problem has no indirect
% reference at all, so G2 is not merely loose -- it is meaningless. Disable it.
floorImposed = isfield(o,'minAltKm') && ~isnan(o.minAltKm);
if floorImposed && ~isempty(tfRef)
    tfRef = [];
    if verbose
        fprintf(['\n  [G2 disabled: an altitude floor is imposed and the indirect\n' ...
                 '   reference solves the UNCONSTRAINED problem, so it is not a\n' ...
                 '   reference for this one.]\n']);
    end
end

%% the measurement
R = dro_residual(o, mu, Tmax, c);
worstKm  = R.RrMax * lStar;                       % POSITION only (see dro_residual)
worstVms = R.RvMax * lStar/p.tStar * 1000;        % velocity, m/s

r2 = vecnorm(o.X(1:3,:) - [1-mu;0;0], 2, 1);
altKm = r2*lStar - rMoonKm;
altAtWorstErrKm = altKm(min(R.kWorstR+1, numel(altKm)));
worstAltKm = altAtWorstErrKm;   % kept for callers; it is the altitude AT the
                                % worst-error interval, not the trajectory minimum

%% the gates
G = struct('id',{},'name',{},'value',{},'tol',{},'pass',{},'units',{});
add = @(id,nm,val,tol,ps,un) struct('id',id,'name',nm,'value',val,'tol',tol,'pass',ps,'units',un);

G(end+1) = add('G1','local POSITION accuracy (worst interval)', worstKm, posTolKm, ...
               worstKm <= posTolKm, 'km');
G(end+1) = add('G1v','local VELOCITY accuracy (worst interval)', worstVms, 1.0, ...
               worstVms <= 1.0, 'm/s');

globKm = NaN;
if ~isempty(globTolKm)
    % shared verifier (migration #4): costate_common/flown_control_error
    if isempty(which('flown_control_error'))
        addpath(fullfile(fileparts(fileparts(fileparts(fileparts( ...
            mfilename('fullpath'))))), 'costate_common'));
    end
    [gr, gv] = flown_control_error(o, mu, Tmax, c);
    globKm  = gr * lStar;
    globVms = gv * lStar/p.tStar * 1000;
    G(end+1) = add('G1b','GLOBAL POSITION accuracy (control flown end to end)', ...
                   globKm, globTolKm, globKm <= globTolKm, 'km');
    G(end+1) = add('G1bv','GLOBAL VELOCITY accuracy (end to end)', globVms, 1.0, ...
                   globVms <= 1.0, 'm/s');
end

if isempty(tfRef)
    tfErrRel = NaN;
else
    tfErrRel = abs(o.tf - tfRef)/tfRef;
    G(end+1) = add('G2','agreement with the indirect reference t_f', tfErrRel, ...
                   tfRelTol, tfErrRel <= tfRelTol, 'rel');
end

G(end+1) = add('G3','NLP defect', o.maxDefect, defectTol, ...
               o.maxDefect <= defectTol, '-');
G(end+1) = add('G4','control unit-norm error', o.maxUnit, 1e-8, ...
               o.maxUnit <= 1e-8, '-');
G(end+1) = add('G5','terminal boundary error', o.termErr, 1e-9, ...
               o.termErr <= 1e-9, '-');
% G6 is ADVISORY. Full throttle is an empirical property of the extremals found
% here, not a theorem: minimum-time problems can admit legitimate coast arcs,
% particularly under path constraints. It is reported, and excluded from passAll.
G(end+1) = add('G6*','throttle saturation, ADVISORY (min u)', ...
               o.thrMin, 1-1e-6, o.thrMin >= 1-1e-6, '-');

if isfield(o,'tfSpread')
    G(end+1) = add('G8','lifted-time spread (all copies equal)', o.tfSpread, ...
                   1e-10, o.tfSpread <= 1e-10, '-');
end
if isfield(o,'maxInterp') && ~isempty(o.maxInterp) && o.maxInterp > 0
    G(end+1) = add('G9','Hermite interpolation residual', o.maxInterp, ...
                   1e-9, o.maxInterp <= 1e-9, '-');
end

% the floor, if one was imposed, checked BETWEEN nodes as well as at them
if isfield(o,'minAltKm') && ~isnan(o.minAltKm)
    trueMinKm = true_min_altitude(o, mu, Tmax, c, lStar, rMoonKm);
    G(end+1) = add('G7','altitude floor honoured between nodes', trueMinKm, ...
                   o.minAltKm, trueMinKm >= o.minAltKm, 'km');
end

% THE VERDICT INCLUDES EVERY MANDATORY GATE. An earlier version keyed the verdict
% on G1 and G2 alone while printing the word CERTIFIED, so a solution could be
% declared certified while failing terminal feasibility or the altitude floor.
% Advisory gates (id ending in '*') are reported but excluded.
isAdv = cellfun(@(x) endsWith(x,'*'), {G.id});
pass    = all([G(~isAdv).pass]);
passAll = all([G.pass]);

C = struct('pass',pass, 'passAll',passAll, 'gates',G, 'resid',R, ...
           'worstAltKm',worstAltKm, 'tfErrRel',tfErrRel, 'worstKm',worstKm, ...
           'globKm',globKm, 'sumLocalKm',sum(R.Rr)*lStar, 'worstVms',worstVms);

%% report
if verbose
    fprintf('\n===== CERTIFICATION: DRO->tulip min-time (%s, N = %d) =====\n', ...
        R.scheme, numel(o.s)-1);
    fprintf('%-4s %-48s %12s %12s  %s\n','','gate','value','tolerance','');
    for k = 1:numel(G)
        fprintf('%-4s %-48s %12.4e %12.4e  %s\n', G(k).id, G(k).name, ...
            G(k).value, G(k).tol, local_mark(G(k).pass));
    end
    fprintf('  worst-error interval sits at %.0f km lunar altitude\n', worstAltKm);
    fprintf('  control reconstruction (%s): min direction norm %.4f, throttle overshoot %.2e\n', ...
        R.scheme, R.minDirNorm, R.thrOvershoot);
    fprintf('  t_f = %.6f ND = %.3f days\n', o.tf, o.tf*p.tStar/86400);
    fprintf('  local POSITION error: max %.4f km, sum %.4f km', worstKm, sum(R.Rr)*lStar);
    if ~isnan(globKm)
        fprintf(', global %.4f km (%s)\n', globKm, ...
            local_accum(globKm, sum(R.Rr)*lStar));
        fprintf('  local VELOCITY error: max %.4f m/s, global %.4f m/s\n', ...
            worstVms, globVms);
    else
        fprintf('\n');
    end
    if pass
        fprintf(['  VERDICT: CERTIFIED -- accurate to the stated tolerances, and\n' ...
                 '           (where checked) agreeing with the reference extremal.\n' ...
                 '           This certifies TRANSCRIPTION ACCURACY and reproduction\n' ...
                 '           of a LOCAL extremal. It does not certify global\n' ...
                 '           optimality, and it says nothing about whether the\n' ...
                 '           underlying OCP is well posed.\n']);
    else
        fprintf('  VERDICT: NOT CERTIFIED -- this t_f is not an answer\n');
    end
end
end

er = norm(z(1:3) - o.X(1:3,end));
ev = norm(z(4:6) - o.X(4:6,end));
end

% ---------------------------------------------------------------------------
function s = local_accum(globKm, sumKm)
% LOCAL_ACCUM  Whether the local errors accumulated or cancelled.
% INPUTS: globKm; sumKm   OUTPUTS: s [char]
if globKm < sumKm, s = 'CANCEL'; else, s = 'ACCUMULATE'; end
end

end



% ---------------------------------------------------------------------------
function s = local_mark(tf)
% LOCAL_MARK  'PASS'/'FAIL' text for the report.
% INPUTS: tf (logical)   OUTPUTS: s [char]
if tf, s = 'PASS'; else, s = 'FAIL'; end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
% INPUTS: s; f; dflt   OUTPUTS: v
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
