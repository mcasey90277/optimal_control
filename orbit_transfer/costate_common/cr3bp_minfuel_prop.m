function [yh, PHI, T, Y] = cr3bp_minfuel_prop(dt, y0, needSTM, Tmax, c, muStar, smooth)
%% Purpose:
%
%   Propagates the smoothed energy->fuel PMP state (and optionally its
%   14x14 STM) for a time dt -- the homotopy sibling of
%   cr3bp_minenergy_prop, in the [yh, PHI] = prop(dt, y0, needSTM) shape
%   ms_bvp expects. Variational equations ride along with A from
%   cr3bp_minfuel_pmp (exact AD). Integrator ode113 at tfMinProp's
%   tolerances (RelTol 1e-10, AbsTol 1e-12).
%
%  ASSUMPTIONS / NOTES:
%
% • Throws (via ode113) on integrator collapse -- ms_bvp converts the
%   throw into a rejected iterate, as required by its contract.
% • The 'huber' family's throttle law is DISCONTINUOUS at Q = 1 (see
%   cr3bp_minfuel_pmp). It is propagated EVENT-SPLIT: each smooth branch
%   with the law pinned, and at every crossing the STM gets the saltation
%   update Phi+ = [I + (F+ - F-) n'/(n'F-)] Phi-, n = grad Q. Without it
%   the shooting Jacobian is wrong on every switching segment (measured
%   2026-09-05: 1.35e-3 vs 1.35e-7 against finite differences;
%   tests/test_huber_saltation). Grazing crossings throw.
%
%% References:
%   [1] Bertrand & Epenoy, "New smoothing techniques for solving bang-bang
%       optimal control problems," OCAM 23(4), 2002 (the 'eps' family).
%   [2] Hairer, Norsett & Wanner, "Solving ODEs I", II.6 (sensitivity
%       across a discontinuity: the saltation / jump matrix).
%
%% Inputs:
%
%  dt                       double                  Propagation time, ND
%
%  y0                       [14 x 1]                Initial PMP state
%
%  needSTM                  logical                 Integrate the STM too
%
%  Tmax, c, muStar          double                  As cr3bp_minfuel_pmp
%
%  smooth                   struct                  {family, p} smoothing
%                                                   spec (cr3bp_minfuel_pmp)
%
%% Outputs:
%
%  yh                       [14 x 1]                State at dt
%
%  PHI                      [14 x 14] or []         STM d yh / d y0
%
%  T                        [nT x 1]                Integrator times
%
%  Y                        [nT x 14]               States at T
%
%% Revision History:
%  M. Casey                                                   (c) 09/02/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: one arc at eps = 0.3:
     mu_ = 0.012150585609624;  T_ = 0.1756418;  c_ = 8.673746;
     y_ = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 0.97; ...
           15.6; 32.9; -0.09; -0.10; 0.045; -0.00015; 0.13];
     [yh_, PHI_] = cr3bp_minfuel_prop(0.4, y_, true, T_, c_, mu_, ...
                                      struct('family','eps','p',0.3));
     fprintf('demo: |r(0.4)| = %.6f, m = %.6f, |PHI| = %.3e\n', ...
             sqrt(sum(yh_(1:3).^2)), yh_(7), max(abs(PHI_(:))));
     return
end

y0 = y0(:);
if dt == 0
    yh = y0;  PHI = eye(14);  T = 0;  Y = y0';
    if ~needSTM, PHI = []; end
    return
end
odeOpts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
if strcmp(smooth.family, 'huber')
    [yh, PHI, T, Y] = propHuber(dt, y0, needSTM, Tmax, c, muStar, smooth, odeOpts);
    return
end
if needSTM
    z0 = [y0; reshape(eye(14), [], 1)];
    [T, Z] = ode113(@(t, z) rhsSTM(z, Tmax, c, muStar, smooth), [0 dt], z0, odeOpts);
    yh  = Z(end, 1:14)';
    PHI = reshape(Z(end, 15:210), 14, 14);
    Y   = Z(:, 1:14);
else
    [T, Y] = ode113(@(t, y) cr3bp_minfuel_pmp(y, Tmax, c, muStar, smooth), ...
                    [0 dt], y0, odeOpts);
    yh  = Y(end, :)';
    PHI = [];
end
end

% ------------------------------------------------------------------------
function [yh, PHI, T, Y] = propHuber(dt, y0, needSTM, Tmax, c, muStar, smooth, odeOpts)
% PROPHUBER  Event-split propagation for the DISCONTINUOUS huber law.
% The throttle jumps at Q = 1 (s: p*Q -> 1). Each smooth branch is
% integrated with the law PINNED (smooth.branch 'lo'|'hi'); at a crossing
% the STM receives the saltation update
%     Phi+ = [ I + (F+ - F-) n' / (n' F-) ] Phi-,   n = grad Q,
% which the branchwise variational equation alone omits (review
% 2026-09-05, P0.2). Grazing crossings (n'F- ~ 0) throw.
% INPUTS: as cr3bp_minfuel_prop + odeOpts.  OUTPUTS: as cr3bp_minfuel_prop.
assert(dt > 0, 'cr3bp_minfuel_prop:backward', ...
       'event-split huber propagation supports dt > 0 only (got %g)', dt);
sLo = smooth;  sLo.branch = 'lo';
sHi = smooth;  sHi.branch = 'hi';
Q0 = switchQ(y0, Tmax, c);
if abs(Q0 - 1) <= 1e-12                          % ON the switch surface: outgoing branch
    hi = cr3bp_minfuel_qdot(y0, Tmax, c) > 0;    % from the (throttle-independent) Qdot
else
    hi = Q0 > 1;
end
t0  = 0;  y = y0;  PHI = eye(14);
T = zeros(0,1);  Y = zeros(0,14);
maxSw = 200;  nSw = 0;
while true
    if hi, sm = sHi;  dirn = -1;  else, sm = sLo;  dirn = +1; end
    eo = odeset(odeOpts, 'Events', @(t, z) evQ(z, Tmax, c, dirn));
    if needSTM
        z0 = [y; reshape(PHI, [], 1)];
        [Ts, Zs, te, ze] = ode113(@(t, z) rhsSTM(z, Tmax, c, muStar, sm), [t0 dt], z0, eo);
        Ys = Zs(:, 1:14);
    else
        [Ts, Ys, te, ze] = ode113(@(t, z) cr3bp_minfuel_pmp(z, Tmax, c, muStar, sm), ...
                                  [t0 dt], y, eo);
        Zs = Ys;
    end
    T = [T; Ts];  Y = [Y; Ys];                     %#ok<AGROW>
    hit = ~isempty(te) && te(end) < dt - 1e-14;
    if ~hit
        yh = Zs(end, 1:14)';
        if needSTM, PHI = reshape(Zs(end, 15:210), 14, 14); else, PHI = []; end
        return
    end
    nSw = nSw + 1;
    assert(nSw <= maxSw, 'cr3bp_minfuel_prop:chatter', ...
           'huber law: more than %d switches in one segment', maxSw);
    t0 = te(end);  ye = ze(end, :)';
    y = ye(1:14);
    if t0 > dt - 1e-9*dt
        warning('cr3bp_minfuel_prop:boundarySwitch', ...
                'huber switch within 1e-9*dt of the segment end (t = %.12g): one-sided STM', t0);
    end
    % Transversality of the crossing = Qdot, which is continuous across the
    % switch (n'F+ == n'F-, measured 9e-16), so it is tested against the
    % size of its own terms, in BOTH propagation modes (Astra review #2):
    Qd = cr3bp_minfuel_qdot(y, Tmax, c);
    qScale = Tmax * sqrt(sum(y(8:10).^2)) / y(7);
    assert(abs(Qd) > 1e-10 * max(qScale, 1e-300), ...
           'cr3bp_minfuel_prop:grazing', ...
           'huber switch is grazing (|Qdot| = %.2e vs scale %.2e)', abs(Qd), qScale);
    if needSTM
        PHI = reshape(ye(15:210), 14, 14);
        Fm = cr3bp_minfuel_pmp(y, Tmax, c, muStar, sm);            % incoming branch
        if hi, smP = sLo; else, smP = sHi; end
        Fp = cr3bp_minfuel_pmp(y, Tmax, c, muStar, smP);           % outgoing branch
        n  = gradQ(y, Tmax, c);
        PHI = (eye(14) + (Fp - Fm) * n' / (n' * Fm)) * PHI;        % n'Fm == Qd
    end
    hi = ~hi;
end
end

function Q = switchQ(y, Tmax, c)
% SWITCHQ  Q = T(|lam_v|/m + lam_m/c) (as cr3bp_minfuel_pmp).
% INPUTS: y [14x1]; Tmax; c.  OUTPUTS: Q double.
rho = sqrt(y(11)^2 + y(12)^2 + y(13)^2 + 1e-300);
Q = Tmax * (rho / y(7) + y(14) / c);
end

function n = gradQ(y, Tmax, c)
% GRADQ  dQ/dy [14x1]: nonzero in m (7), lam_v (11:13), lam_m (14).
% INPUTS: y [14x1]; Tmax; c.  OUTPUTS: n [14x1].
rho = sqrt(y(11)^2 + y(12)^2 + y(13)^2 + 1e-300);
n = zeros(14, 1);
n(7)     = -Tmax * rho / y(7)^2;
n(11:13) =  Tmax * y(11:13) / (y(7) * rho);
n(14)    =  Tmax / c;
end

function [val, isterm, dirn] = evQ(z, Tmax, c, dirn)
% EVQ  Terminal event Q - 1 = 0, crossing in direction dirn only.
% INPUTS: z [14 or 210 x1]; Tmax; c; dirn.  OUTPUTS: val; isterm; dirn.
val = switchQ(z(1:14), Tmax, c) - 1;
isterm = 1;
end

% ------------------------------------------------------------------------
function dz = rhsSTM(z, Tmax, c, muStar, smooth)
% RHSSTM  PMP state + STM variational RHS.  INPUTS: z [210x1]; field
% params + smooth.  OUTPUTS: dz [210x1].
[F, A] = cr3bp_minfuel_pmp(z(1:14), Tmax, c, muStar, smooth);
PHI = reshape(z(15:210), 14, 14);
dz = [F; reshape(A*PHI, [], 1)];
end
