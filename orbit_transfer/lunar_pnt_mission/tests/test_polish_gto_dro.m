function ok = test_polish_gto_dro()
%% Purpose:
%
%   Verifies polish_gto_dro (leg 1, GTO -> DRO) on the clean perigee cell
%   (iD 1, iA 1) of the first sheet: (1) identity -- the sheet orientation
%   and the cell's own arrival phase return the catalog t_f (1e-8) and z8
%   (1e-6); (2) +10 deg and a slot 0.03 of a period past the cell's
%   arrival phase -- accepted, flown miss < 0.5 km and < 1 cm/s, departure
%   state exactly L.gtoState(orientDeg, sD), arrival phase
%   mod(slotPhaseAtDep + tf/tauD, 1), the arrival state checked against an
%   INDEPENDENT tight propagation of the DRO from its phase-0 state, and
%   the SLIDE CONVENTION checked independently (the slot coasted
%   ballistically from slotPhaseAtDep for P.tf on this test's own CR3BP
%   right-hand side lands on P.rvf, 1e-6 ND); familyChanged equals
%   |tf - tfCell|/tfCell > familyFrac;
%   (3) Review Focus 1 -- a 60 deg target in ONE step with no halvings
%   (measured 2026-11-11: ms_tfmin hits 100 iterations at residual 1.5)
%   returns ok = false, a reason, orientReached at the sheet orientation
%   and an empty z8; (3b) deliberate branch-jump failure (jumpFrac 1e-4,
%   the 2 deg step moves t_f by 1.3e-3 of itself); (4) deliberate failure
%   of the miss gate through opts.gateOverride, and an unknown override
%   name refused.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  final fix wave: miss < 0.5 km, independent slide
%            convention, familyChanged                        (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
assert(strcmp(fileparts(which('polish_gto_dro')), fullfile(fileparts(here), 'lib')), ...
       'polish_gto_dro not from this campaign');

L  = leg1_inputs(struct());
iS = 1;  iD = 1;  iA = 1;
assert(L.clean(iD, iA, iS), 'cell (%d,%d) of sheet %d is not clean', iD, iA, iS);
o0     = L.sheetOrient(iS);
tfCell = L.tf(iD, iA, iS);
z8Cell = L.z8(:, iD, iA, iS);
sD     = L.sD(iD);
base   = struct('iD', iD, 'iA', iA, 'iS', iS, 'tDep', 0.25);

% ---- 1. identity ----
c1 = base;
c1.orientDeg      = o0;
c1.slotPhaseAtDep = mod(L.sA(iA) - tfCell/L.tauD, 1);
P1 = polish_gto_dro(L, c1, struct());
fprintf(['  identity: ok %d (%s)  tf %.12f vs cell %.12f (%.1e)  |dz8| %.1e  ' ...
         'nEval %d  miss %.3g km %.3g m/s\n'], P1.ok, P1.reason, P1.tf, tfCell, ...
        abs(P1.tf - tfCell), max(abs(P1.z8 - z8Cell)), P1.nEval, P1.missKm, P1.missMps);
assert(P1.ok, 'identity not ok: %s', P1.reason);
assert(abs(P1.tf - tfCell) < 1e-8, 'identity tf off by %.3g', abs(P1.tf - tfCell));
assert(max(abs(P1.z8 - z8Cell)) < 1e-6, 'identity z8 off by %.3g', max(abs(P1.z8 - z8Cell)));
assert(isequal(P1.walk.orientDeg, o0) && P1.orientReached == o0, 'identity walked');

% ---- 2. +10 deg, slot 0.03 past the arrival phase ----
c2 = base;
c2.orientDeg      = o0 + 10;
c2.slotPhaseAtDep = mod(L.sA(iA) + 0.03 - tfCell/L.tauD, 1);
P2 = polish_gto_dro(L, c2, struct('verbose', true));
fprintf(['  +10 deg / +0.03: ok %d (%s)  tf %.9f (cell %.9f)  arrPhase %.6f  ' ...
         'miss %.3g km %.3g m/s  res %.2e  nEval %d  walk %d steps  %.1f s\n'], ...
        P2.ok, P2.reason, P2.tf, tfCell, P2.arrPhase, P2.missKm, P2.missMps, P2.res, ...
        P2.nEval, numel(P2.walk.orientDeg) - 1, P2.wallS);
assert(P2.ok, '+10 deg not ok: %s', P2.reason);
assert(P2.missKm < 0.5 && P2.missMps < 0.01, 'flown miss %.3g km %.3g m/s', P2.missKm, P2.missMps);
assert(P2.res < 1e-8, 'residual %.3g', P2.res);
assert(max(abs(P2.rv0 - L.gtoState(c2.orientDeg, sD))) < 1e-12, 'departure state');
assert(P2.orientDeg == c2.orientDeg && P2.orientReached == c2.orientDeg && ...
       P2.walk.orientDeg(end) == c2.orientDeg && P2.walk.orientDeg(1) == o0, 'walk ends');
assert(abs(P2.arrPhase - mod(c2.slotPhaseAtDep + P2.tf/L.tauD, 1)) < 1e-12, 'arrival phase');
assert(max(abs(P2.rvf - L.slotState(P2.arrPhase))) < 1e-12, 'arrival state vs slotState');
assert(abs(P2.tArr - (c2.tDep + P2.tf)) < 1e-12 && P2.tDep == c2.tDep, 'epochs');
assert(abs(P2.sD - sD) < 1e-15, 'sD');
% INDEPENDENT arrival check: the DRO flown by ode113 on pumpkyn's CR3BP
% field from its phase-0 state for arrPhase*tauD
tProp = P2.arrPhase*L.tauD;
[~, Yd] = ode113(@pumpkyn.cr3bp.eom, [0 0.5*tProp tProp], L.slotState(0).', ...
                 odeset('RelTol', 1e-13, 'AbsTol', 1e-13), L.muStar, 1);
dArr = max(abs(P2.rvf - Yd(end, :).'));
fprintf('  arrival state vs direct DRO propagation (%.4f ND): %.2e ND (%.3f km)\n', ...
        tProp, dArr, dArr*L.lStar);
assert(dArr < 1e-7, 'arrival state differs from the propagated DRO by %.3g ND', dArr);
% INDEPENDENT SLIDE CONVENTION: the slot, at DRO phase slotPhaseAtDep at the
% departure, coasted ballistically for P.tf ND on this test's OWN CR3BP
% right-hand side, is where the polished transfer arrives.
% Measured 2026-11-11: started from L.slotState(slotPhaseAtDep) the coast
% lands 1.41e-6 ND (0.55 km) from rvf -- the floor of the DRO phase_state
% handle, which is 9.2e-8 ND off the true orbit at phase 0.4846 (105
% stored points; 3.8e-6 ND after one period from phase 0.5, 2e-12 from
% phase 0) and that start error grows over 0.545 ND. That number is
% PRINTED. The ASSERTION (1e-6 ND) starts the slot where the handle is
% exact, at phase 0, and coasts it slotPhaseAtDep*tauD to the departure
% (the slot's departure-epoch state, by propagation) and then P.tf more:
% every phase advance is done by the ODE, none by the library's phase
% arithmetic. A wrong convention (slot pinned at its departure phase, or
% the phase counted from the arrival epoch) misses by 0.03 to 0.5 of a
% period.
op13 = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[~, Ys] = ode113(@(t, x) cr3bp_rhs(x, L.muStar), [0 0.5*P2.tf P2.tf], L.slotState(c2.slotPhaseAtDep), op13);
dSlideH = max(abs(P2.rvf - Ys(end, :).'));
tDepSlot = c2.slotPhaseAtDep*L.tauD;
[~, Yc] = ode113(@(t, x) cr3bp_rhs(x, L.muStar), [0 tDepSlot tDepSlot + P2.tf], L.slotState(0), op13);
dSlide = max(abs(P2.rvf - Yc(end, :).'));
fprintf(['  slide convention: slot coasted %.6f ND from its departure phase %.6f lands %.2e ND (%.3f km) ' ...
         'from rvf when started from the handle (its interpolation floor), %.2e ND (%.2g km) when the ' ...
         'departure state is itself propagated from phase 0\n'], P2.tf, c2.slotPhaseAtDep, dSlideH, ...
        dSlideH*L.lStar, dSlide, dSlide*L.lStar);
assert(dSlide < 1e-6, 'slide convention: the coasted slot misses rvf by %.3g ND', dSlide);
% family flag: |tf - tfCell|/tfCell > familyFrac (0.02 by default)
assert(islogical(P2.familyChanged) && P2.familyChanged == (abs(P2.tf - P2.tfCell)/P2.tfCell > 0.02) && ...
       ~P1.familyChanged, 'familyChanged: P1 %d, P2 %d (rel %.3g)', P1.familyChanged, P2.familyChanged, ...
       abs(P2.tf - P2.tfCell)/P2.tfCell);
fprintf('  familyChanged: identity %d, +10 deg %d (|tf - tfCell|/tfCell = %.3g)\n', P1.familyChanged, ...
        P2.familyChanged, abs(P2.tf - P2.tfCell)/P2.tfCell);

% ---- 3. Review Focus 1: the walk loses the solution ----
c3 = base;
c3.orientDeg      = o0 + 60;
c3.slotPhaseAtDep = c1.slotPhaseAtDep;
P3 = polish_gto_dro(L, c3, struct('stepDeg', 60, 'maxHalvings', 0));
fprintf('  60 deg in one step: ok %d (%s)  orientReached %g\n', P3.ok, P3.reason, P3.orientReached);
assert(~P3.ok && ~isempty(P3.reason) && contains(P3.reason, 'orientation walk stalled'), ...
       'stall: %s', P3.reason);
assert(P3.orientReached == o0, 'stall: orientReached %g', P3.orientReached);
assert(isempty(P3.z8) && isnan(P3.tf), 'stall: a transfer was returned');

% ---- 3b. deliberate branch-jump failure ----
P3b = polish_gto_dro(L, c2, struct('jumpFrac', 1e-4, 'maxHalvings', 0));
fprintf('  jumpFrac 1e-4: ok %d (%s)  orientReached %g\n', P3b.ok, P3b.reason, P3b.orientReached);
assert(~P3b.ok && contains(P3b.reason, 'branch jump') && P3b.orientReached == o0 && ...
       isempty(P3b.z8), 'branch jump: %s', P3b.reason);

% ---- 4. deliberate failure of the miss gate ----
P4 = polish_gto_dro(L, c1, struct('gateOverride', struct('flyKmMax', 1e-9)));
fprintf('  flyKmMax 1e-9: ok %d (%s)  miss %.3g km\n', P4.ok, P4.reason, P4.missKm);
assert(~P4.ok && startsWith(P4.reason, 'gate failed') && contains(P4.reason, 'fly'), ...
       'miss gate: %s', P4.reason);
assert(isfinite(P4.missKm) && isempty(P4.z8), 'miss gate: shape');
assert(isequal(sort(fieldnames(P4)), sort(fieldnames(P1))), 'struct shape differs');
threw = false;
try
    polish_gto_dro(L, c1, struct('gateOverride', struct('noSuchGate', 1)));
catch err
    threw = strcmp(err.identifier, 'polish_gto_dro:gateOverride');
end
assert(threw, 'unknown override name not refused');

fprintf('  PASS\n');
ok = true;
end

% ------------------------------------------------------------------------
function dx = cr3bp_rhs(x, mu)
% CR3BP_RHS  Ballistic CR3BP, barycentric rotating frame (Earth at -mu, Moon
%   at 1 - mu), written here independently of pumpkyn and the library.
%   INPUTS: x [6 x 1] ND; mu scalar.  OUTPUTS: dx [6 x 1].
r1 = sqrt((x(1) + mu)^2 + x(2)^2 + x(3)^2);
r2 = sqrt((x(1) - 1 + mu)^2 + x(2)^2 + x(3)^2);
dx = [x(4:6);
      2*x(5) + x(1) - (1 - mu)*(x(1) + mu)/r1^3 - mu*(x(1) - 1 + mu)/r2^3;
      -2*x(4) + x(2) - (1 - mu)*x(2)/r1^3 - mu*x(2)/r2^3;
      -(1 - mu)*x(3)/r1^3 - mu*x(3)/r2^3];
end
