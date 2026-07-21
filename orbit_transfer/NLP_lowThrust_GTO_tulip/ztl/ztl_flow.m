function out = ztl_flow(y0, tspan, P, wantSTM)
% ZTL_FLOW  Integrate the ramp-family PMP flow (+ variational STM) through
% the 3-regime automaton: segment-by-segment between regime-boundary events.
%
% Regimes of the BE ramp throttle (boundaries S = -eps and S = +eps):
%   ON (S<=-eps, u=1) | MEDIUM (|S|<eps, ramp) | OFF (S>=+eps, u=0)
% Integration NEVER steps across a boundary: each segment runs ode89 with a
% DIRECTIONAL terminal event on the boundary being approached, then restarts
% in the new regime (directional events cannot re-fire at the restart point).
% For eps>0 the field is continuous across boundaries (no STM correction);
% for eps=0 the field JUMPS at S=0 and the saltation matrix
%   Psi = I + (f_plus - f_minus) * dSdy / Sdot_minus
% is composed into the STM. |Sdot| below P.grazeFloor at an event flags a
% graze (Psi untrustworthy).
%
% The STM solves the variational ODE Phi_dot = A(y(t))*Phi with A from
% ztl_A (complex step of the FIELD -- exact; see Z0_BUILD.md par.3).
%
% INPUTS:
%   y0      - initial augmented state [14x1]
%   tspan   - [t0 tf] (ND)
%   P       - struct: .muStar .c .Tmax .eps, optional .odeRelTol [1e-13]
%             .odeAbsTol [1e-15] .grazeFloor [1e-4] .maxSegs [400]
%   wantSTM - (optional) propagate the 14x14 STM [default false]
%
% OUTPUTS:
%   out - struct:
%     .yf     final state [14x1]        .PHI    STM(tf,t0) [14x14] (I if off)
%     .events struct array: t, S, Sdot, from, to, grazed
%     .nSegs  segments used             .flag   0 ok | 1 graze | 2 maxSegs
%     .t,.y   concatenated solver output (diagnostics; row-per-step)
%
% REFERENCES:
%   [1] Zhang et al., JGCD 38(8), 2015 (switch-detected integration + STM).
%   [2] Z0_BUILD.md par.6 (this function's spec).

if nargin < 4, wantSTM = false; end
relTol = getdef(P, 'odeRelTol', 1e-13);
absTol = getdef(P, 'odeAbsTol', 1e-15);
grazeFloor = getdef(P, 'grazeFloor', 1e-4);
maxSegs    = getdef(P, 'maxSegs', 400);
eps_ = P.eps;

t  = tspan(1);  tf = tspan(2);
y  = y0(:);
PHI = eye(14);
events = struct('t',{},'S',{},'Sdot',{},'from',{},'to',{},'grazed',{},'yEv',{});
tAll = [];  yAll = [];
flag = 0;  nSegs = 0;

% initial regime (boundary ties broken by Sdot direction)
regime = classify_regime(y, P);

while t < tf - 1e-13
    nSegs = nSegs + 1;
    if nSegs > maxSegs, flag = 2; break; end

    if wantSTM
        z0 = [y; PHI(:)];
        rhs = @(tt, z) rhs_stm(z, P, regime);
    else
        z0 = y;
        rhs = @(tt, z) ztl_eom(z, P, regime);
    end
    opts = odeset('RelTol', relTol, 'AbsTol', absTol, ...
                  'Events', @(tt, z) boundary_events(z, P, regime));
    [T, Z, ~, ~, IE] = ode89(rhs, [t tf], z0, opts);

    t = T(end);
    zEnd = Z(end, :).';
    y = zEnd(1:14);
    if wantSTM, PHI = reshape(zEnd(15:end), 14, 14); end
    tAll = [tAll; T];  yAll = [yAll; Z(:, 1:14)]; %#ok<AGROW>

    if t >= tf - 1e-13, break; end            % reached final time
    assert(~isempty(IE), 'ztl_flow: segment ended before tf without an event');

    % --- regime transition at the event -------------------------------------
    [fMinus, auxM] = ztl_eom(y, P, regime);
    newRegime = next_regime(regime, IE(end), eps_);
    grazed = abs(auxM.Sdot) < grazeFloor;
    if grazed, flag = max(flag, 1); end

    if eps_ == 0 && wantSTM
        fPlus = ztl_eom(y, P, newRegime);
        dSdy  = dS_dy(y, P);
        Psi   = eye(14) + ((fPlus - fMinus) * dSdy) / auxM.Sdot;
        PHI   = Psi * PHI;
    end
    events(end+1) = struct('t', t, 'S', real(auxM.S), 'Sdot', real(auxM.Sdot), ...
        'from', regime, 'to', newRegime, 'grazed', grazed, 'yEv', y); %#ok<AGROW>
    regime = newRegime;
end

out = struct('yf', y, 'PHI', PHI, 'events', events, 'nSegs', nSegs, ...
             'flag', flag, 't', tAll, 'y', yAll);
end

% ---------------------------------------------------------------------------
function dz = rhs_stm(z, P, regime)
% RHS for [y; Phi(:)]: field + variational equations Phi_dot = A(y)*Phi.
y   = z(1:14);
Phi = reshape(z(15:end), 14, 14);
f   = ztl_eom(y, P, regime);
A   = ztl_A(y, P, regime);
dz  = [f; reshape(A*Phi, [], 1)];
end

function [value, isterminal, direction] = boundary_events(z, P, regime)
% Directional terminal events on the CURRENT regime's exit boundaries.
m = z(7);  lam_v = z(11:13);  lam_m = z(14);
S = 1 - sqrt(sum(lam_v.^2))*P.c/m - lam_m;
eps_ = P.eps;
switch regime
    case 'on'                       % exit upward through S = -eps
        value = S + eps_;  direction = +1;  isterminal = 1;
    case 'off'                      % exit downward through S = +eps
        value = S - eps_;  direction = -1;  isterminal = 1;
    case 'medium'                   % exit either boundary
        value = [S - eps_; S + eps_];
        direction = [+1; -1];  isterminal = [1; 1];
end
end

function regime = classify_regime(y, P)
% Regime from S, boundary ties broken by the sign of Sdot (entering side).
[~, aux] = ztl_eom(y, P, 'off');    % S/Sdot are regime-independent... (Sdot
% depends on u through mDot/lvDot? lvDot no; mDot & lmDot yes -> compute
% Sdot consistently AFTER picking a candidate regime below when on a tie.)
S = real(aux.S);  eps_ = P.eps;
tol = 1e-12;
if S <= -eps_ - tol
    regime = 'on';
elseif S >= eps_ + tol
    regime = 'off';
elseif eps_ > 0 && abs(S) < eps_ - tol
    regime = 'medium';
else
    % on a boundary: pick the side S is moving into (u continuous at eps>0,
    % so Sdot from either neighboring regime agrees to O(tol); at eps=0 use
    % the ON-side Sdot, matching the saltation convention)
    if eps_ == 0
        [~, auxOn] = ztl_eom(y, P, 'on');
        if real(auxOn.Sdot) > 0, regime = 'off'; else, regime = 'on'; end
    else
        candidate = 'medium';
        [~, auxC] = ztl_eom(y, P, candidate);
        sd = real(auxC.Sdot);
        if abs(S - eps_) < abs(S + eps_)   % at the +eps boundary
            if sd > 0, regime = 'off'; else, regime = 'medium'; end
        else                               % at the -eps boundary
            if sd < 0, regime = 'on'; else, regime = 'medium'; end
        end
    end
end
end

function newRegime = next_regime(regime, ie, eps_)
% Regime entered when the ie-th event of the CURRENT regime fires.
switch regime
    case 'on',  newRegime = tern(eps_ > 0, 'medium', 'off');
    case 'off', newRegime = tern(eps_ > 0, 'medium', 'on');
    case 'medium'
        if ie == 1, newRegime = 'off'; else, newRegime = 'on'; end
end
end

function row = dS_dy(y, P)
% dS/dy [1x14] for S = 1 - ||lam_v|| c/m - lam_m (analytic).
m = y(7);  lam_v = y(11:13);
lamvMag = sqrt(sum(lam_v.^2));
row = zeros(1, 14);
row(7)      =  P.c*lamvMag/m^2;
row(11:13)  = -(P.c/m)*(lam_v.'/lamvMag);
row(14)     = -1;
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end
