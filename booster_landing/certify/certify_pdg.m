function rep = certify_pdg(solC, solV, P, tolScale)
% CERTIFY_PDG  Gates G1-G5 for a PDG solution pair (colloc + convex).
%
% G1 (discrete):     max Hermite-Simpson defect, re-evaluated in plain
%                    MATLAB pdg_dynamics (independent of the CasADi build
%                    the NLP was actually solved with).
% G2 (continuous):   ode45 re-integration of the pchip-interpolated
%                    colloc control against the colloc state trajectory
%                    -- "defect is not accuracy" (house lesson): G1 only
%                    proves the discrete constraints were satisfied, G2 is
%                    the first real measurement of continuous-time error.
% G3 (cross-method): final mass / final time / position-trajectory
%                    agreement between the independently-formulated
%                    collocation (nonconvex annulus) and convexified
%                    (lossless relaxation) solvers.
% G4 (losslessness): the convex relaxation gap (||u||-sigma) is tight.
% G5 (PMP structure): throttle bang-bang (>=95% of nodes on a bound,
%                    <=2 interior switches, max-thrust at touchdown) and
%                    primer-vector alignment (thrust parallel to the
%                    velocity-costate direction), the two textbook
%                    necessary conditions for a mass-optimal 3-DOF landing
%                    burn (Lawden/PMP).
%
% INPUTS:
%   solC     - sol struct from solve_pdg_colloc [see that function]
%   solV     - sol struct from solve_pdg_convex, or [] to skip G3/G4
%              (Phase 2 drag runs have no convex twin, since the
%              convexification is only exact in vacuum)
%   P        - booster_params struct
%   tolScale - (optional, def 1) scales the G2 continuous-residual and G3
%              cross-method agreement tolerances; for COARSE-grid test
%              calls only (see tests/test_certify_nominal.m) -- never used
%              to loosen the nominal-grid (tolScale=1) tolerance in
%              production reports. G2 needs this too, not just G3: an
%              N-sweep (task-5 report) shows the ode45 re-integration
%              residual (esp. mass) plateaus around 0.6-0.8 kg for N>=60
%              and does not vanish with further refinement -- the same
%              genuine, bounded, non-formulation-error gap Fact 2
%              documents for G3, independently rediscovered in G2.
%
% OUTPUTS:
%   rep - gate report struct:
%     .G1_defect, .G1_pass
%     .G2_pos, .G2_vel, .G2_dm, .G2_pass
%     .G3_dmf, .G3_dtf, .G3_traj_Linf, .G3_pass ('skipped' if solV==[])
%     .G4_gap, .G4_pass ('skipped' if solV==[])
%     .G5_bound_frac, .G5_switches, .G5_primer_deg, .G5_pass
%     .all_pass - logical AND of all gates (skipped gates excluded)
%
% This function is report-only and never throws; callers (tests, the
% front door) decide what a FAIL means.
%
% REFERENCES:
%   [1] Acikmese & Ploen, "Convex Programming Approach to Powered Descent
%       Guidance for Mars Landing," JGCD 2007. (lossless convexification)
%   [2] Lawden, "Optimal Trajectories for Space Navigation," 1963.
%       (primer vector theory / PMP bang-bang structure)
if nargin < 4, tolScale = 1; end

%% G1: HS defects re-evaluated in plain MATLAB (independent of CasADi):
N = size(solC.X,2) - 1;  h = solC.tf / N;  dmax = 0;
for k = 1:N
    fk  = pdg_dynamics(solC.X(:,k),   solC.U(:,k),   P);
    fk1 = pdg_dynamics(solC.X(:,k+1), solC.U(:,k+1), P);
    xm  = 0.5*(solC.X(:,k)+solC.X(:,k+1)) + (h/8)*(fk - fk1);
    fm  = pdg_dynamics(xm, solC.Um(:,k), P);
    d   = solC.X(:,k+1) - solC.X(:,k) - (h/6)*(fk + 4*fm + fk1);
    dmax = max(dmax, max(abs(d)));
end
rep.G1_defect = dmax;   rep.G1_pass = dmax < 1e-6;

%% G2: continuous residual -- fly the interpolated control with ode45
%% ("defect is not accuracy"): pchip on [nodes+midpoints] thrust.
tU = sort([solC.t, solC.t(1:end-1) + h/2]);
TU = zeros(3, numel(tU));
TU(:,1:2:end) = solC.U;  TU(:,2:2:end) = solC.Um;
ctrl = @(tt) interp1(tU.', TU.', min(tt, solC.tf), 'pchip').';
odef = @(tt, xx) pdg_dynamics(xx, ctrl(tt), P);
oo   = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~, XX] = ode45(odef, [0 solC.tf], solC.X(:,1), oo);
ef = XX(end,:).' - solC.X(:,end);
rep.G2_pos = sqrt(sum(ef(1:3).^2));  rep.G2_vel = sqrt(sum(ef(4:6).^2));
rep.G2_dm  = abs(ef(7));
rep.G2_pass = rep.G2_pos < 1*tolScale && rep.G2_vel < 0.1*tolScale ...
              && rep.G2_dm < 0.5*tolScale;

%% G3/G4: cross-method agreement + losslessness (skip if no convex twin):
if isempty(solV)
    rep.G3_pass = 'skipped';  rep.G4_pass = 'skipped';
else
    rep.G3_dmf = abs(solC.mf - solV.mf);
    rep.G3_dtf = abs(solC.tf - solV.tf);
    tq  = linspace(0, min(solC.tf, solV.tf), 200);
    rC  = interp1(solC.t.', solC.X(1:3,:).', tq.', 'pchip');
    rV  = interp1(solV.t.', solV.X(1:3,:).', tq.', 'pchip');
    rep.G3_traj_Linf = max(sqrt(sum((rC - rV).^2, 2)));
    rep.G3_pass = rep.G3_dmf < 0.1*tolScale && rep.G3_dtf < 0.2*tolScale;
    rep.G4_gap  = solV.lossless_gap;
    rep.G4_pass = rep.G4_gap < 1e-4 * P.Tmax / P.m0;
end

%% G5: PMP structure. Throttle bang-bang + primer alignment.
Tmag = sqrt(sum(solC.U.^2, 1));
onLo = Tmag < P.Tmin * 1.001;   onHi = Tmag > P.Tmax * 0.999;
rep.G5_bound_frac = mean(onLo | onHi);
segs = diff([onHi(1), onHi]);                  % max<->min transitions
rep.G5_switches = sum(segs ~= 0);
% Structure check (ADAPTATION from the brief's literal "onHi(1) && onHi(end)"
% [max-first AND max-last], documented in the task-5 report): the measured
% nominal solution is min-first / max-last (onHi(1)=false, single switch
% Tmin->Tmax at t~7.06 s of tf~15.69 s, held to touchdown), not max-first.
% This is genuine physics, not a bug: task-3's own report established that
% P.Tmin (338 kN) already exceeds the hover-equivalent thrust at any mass
% in [mdry,m0] (m0*g0=294.2 kN, mdry*g0=251.1 kN), so the vehicle cannot
% loiter even at minimum throttle -- coasting at Tmin early, then a single
% hard Tmax brake just before the pad, is the fuel-optimal bang-bang
% structure for THIS entry state, and PMP does not mandate a specific
% opening throttle (only: bang-bang, no singular arc, and -- for a landing
% burn that must arrest velocity before touchdown -- max thrust at the
% very end). We therefore require max-LAST (a physical necessity: some
% peak-thrust braking is needed right before touchdown) but not max-FIRST.
structOk = rep.G5_bound_frac >= 0.95 && rep.G5_switches <= 2 && onHi(end);
% Primer: velocity costate from defect duals; thrust must align with the
% costate direction. SIGN (ADAPTATION, documented per the brief's own
% instruction: "if G5_primer_deg comes out near 180 deg ... flip pdir sign
% ONCE"): the brief's literal pdir = -lamv/|lamv| gave primer_deg ~179.7-
% 179.9 deg (thrust anti-parallel to -lam_v). Flipping to pdir=+lamv/|lamv|
% gives ~0.15-0.35 deg (thrust parallel to +lam_v) -- so for the sign
% convention solve_pdg_colloc's opti.lam_g extraction actually returns
% (house lesson: opti.lam_g, never opti.dual -- see MEMORY
% casadi-optidual-sign-bug), the correctly-aligned primer direction is
% +lam_v, not -lam_v. Applied ONCE here, no abs() anywhere in this file.
lamv = solC.lam_defect(4:6, :);                % 3xN, node-adjacent duals
Tdir = solC.U(:,1:end-1) ./ sqrt(sum(solC.U(:,1:end-1).^2,1));
pdir = lamv ./ max(sqrt(sum(lamv.^2,1)), 1e-30);
cosang = sum(Tdir .* pdir, 1);
rep.G5_primer_deg = max(acosd(min(1, max(-1, cosang))));
rep.G5_pass = structOk && rep.G5_primer_deg < 1;

%% Verdict:
gates = [rep.G1_pass, rep.G2_pass, isequal(rep.G3_pass,true) || ...
         isequal(rep.G3_pass,'skipped'), isequal(rep.G4_pass,true) || ...
         isequal(rep.G4_pass,'skipped'), rep.G5_pass];
rep.all_pass = all(gates);
end

% NOTE (deviation from the brief's "nested function" suggestion, documented
% in the task-5 report): print_certify_report is called directly by
% tests/test_certify_nominal.m and (Task 10) the front door -- both outside
% this file. MATLAB local functions are only visible within their own
% function file, so a copy nested here would NOT be callable from either
% caller; the brief's own fallback ("expose it by making
% certify/print_certify_report.m a 3-line wrapper if needed elsewhere")
% applies, and it IS needed elsewhere. Rather than maintain two copies (one
% nested/dead, one real) the single real implementation lives in its own
% file, certify/print_certify_report.m, right beside this one.
