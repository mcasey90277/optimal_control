function out = ztl_flow_sun(Y0, tauF, sigSpan, P, wantSTM)
% ZTL_FLOW_SUN  Integrate the Sundman-regularized PMP flow (+ STM + tauF
% sensitivity) over a sigma-interval, through the 3-regime automaton.
%
% Independent variable sigma; dY/dsigma = tauF * G(Y), G = ztl_eom_sun RHS
% (dY/dtau). When wantSTM, also propagates Phi = dY/dY0 (15x15) and
% w = dY/dtauF (15x1, INHOMOGENEOUS: dw/dsigma = tauF*A*w + G). Regime events
% on S(y) as in ztl_flow (directional; saltation at eps=0).
%
% INPUTS:
%   Y0      - initial Sundman state [15x1]: [r;v;m;lam_r;lam_v;lam_m;t]
%   tauF    - total Sundman length (scalar; scales the RHS)
%   sigSpan - [sigma0 sigma1]
%   P       - .muStar .c .Tmax .eps .pSund, opt .odeRelTol[1e-13]
%             .odeAbsTol[1e-15] .grazeFloor[1e-4] .maxSegs[400]
%   wantSTM - propagate Phi and w [default false]
%
% OUTPUTS:
%   out - .Yf [15x1] .PHI [15x15] .w [15x1] (dYf/dtauF; 0 if ~wantSTM)
%         .events (t,S,Sdot,from,to,grazed) .nSegs .flag (0|1 graze|2 maxSegs)
%         .sig,.Y (diagnostics)
%
% REFERENCES: ztl_flow.m (physical-time analog); SUN_BUILD.md.

if nargin < 5, wantSTM = false; end
relTol = getdef(P,'odeRelTol',1e-13);  absTol = getdef(P,'odeAbsTol',1e-15);
grazeFloor = getdef(P,'grazeFloor',1e-4);  maxSegs = getdef(P,'maxSegs',400);
eps_ = P.eps;

s = sigSpan(1);  sf = sigSpan(2);
Y = Y0(:);  PHI = eye(15);  w = zeros(15,1);
events = struct('t',{},'S',{},'Sdot',{},'from',{},'to',{},'grazed',{},'Yev',{});
sAll = [];  YAll = [];  flag = 0;  nSegs = 0;
regime = classify_regime(Y, P);

while s < sf - 1e-13
    nSegs = nSegs + 1;
    if nSegs > maxSegs, flag = 2; break; end
    if wantSTM
        z0 = [Y; PHI(:); w];
        rhs = @(ss,z) rhs_stm_sun(z, P, regime, tauF);
    else
        z0 = Y;
        rhs = @(ss,z) tauF*ztl_eom_sun(z, P, regime);
    end
    opts = odeset('RelTol',relTol,'AbsTol',absTol, ...
                  'Events',@(ss,z) boundary_events(z, P, regime));
    [Sg, Z, ~, ~, IE] = ode89(rhs, [s sf], z0, opts);
    s = Sg(end);  zEnd = Z(end,:).';
    Y = zEnd(1:15);
    if wantSTM, PHI = reshape(zEnd(16:240),15,15);  w = zEnd(241:255); end
    sAll = [sAll; Sg];  YAll = [YAll; Z(:,1:15)]; %#ok<AGROW>
    if s >= sf - 1e-13, break; end
    assert(~isempty(IE), 'ztl_flow_sun: segment ended before sigma_f without an event');

    [Gminus, auxM] = ztl_eom_sun(Y, P, regime);
    newRegime = next_regime(regime, IE(end), eps_);
    SdotTau = auxM_Sdot(Y, P, Gminus);
    grazed = abs(SdotTau) < grazeFloor;
    if grazed, flag = max(flag,1); end
    if eps_ == 0 && wantSTM
        Gplus = ztl_eom_sun(Y, P, newRegime);
        dSdY = dS_dY(Y, P);
        Psi = eye(15) + ((Gplus - Gminus) * dSdY) / SdotTau;
        PHI = Psi * PHI;  w = Psi * w;
    end
    events(end+1) = struct('t',real(Y(15)),'S',real(auxM.S),'Sdot',real(SdotTau), ...
        'from',regime,'to',newRegime,'grazed',grazed,'Yev',Y); %#ok<AGROW>
    regime = newRegime;
end

out = struct('Yf',Y,'PHI',PHI,'w',w,'events',events,'nSegs',nSegs, ...
             'flag',flag,'sig',sAll,'Y',YAll);
end

% ---------------------------------------------------------------------------
function dz = rhs_stm_sun(z, P, regime, tauF)
Y = z(1:15);  Phi = reshape(z(16:240),15,15);  w = z(241:255);
G = ztl_eom_sun(Y, P, regime);
A = ztl_A_sun(Y, P, regime);
dz = [tauF*G; reshape(tauF*(A*Phi),[],1); tauF*(A*w) + G];
end

function [value,isterminal,direction] = boundary_events(z, P, regime)
m = z(7);  lam_v = z(11:13);  lam_m = z(14);
S = 1 - sqrt(sum(lam_v.^2))*P.c/m - lam_m;
eps_ = P.eps;
switch regime
    case 'on',     value = S + eps_;  direction = +1;  isterminal = 1;
    case 'off',    value = S - eps_;  direction = -1;  isterminal = 1;
    case 'medium', value = [S - eps_; S + eps_];  direction = [+1;-1];  isterminal = [1;1];
end
end

function SdotTau = auxM_Sdot(Y, P, G)
% dS/dtau = dS/dY * (dY/dtau) = dS_dY * G   (S depends on y only).
SdotTau = dS_dY(Y, P) * G;
end

function regime = classify_regime(Y, P)
m = Y(7);  lam_v = Y(11:13);  lam_m = Y(14);  eps_ = P.eps;  tol = 1e-12;
S = 1 - sqrt(sum(lam_v.^2))*P.c/m - lam_m;
if S <= -eps_ - tol
    regime = 'on';
elseif S >= eps_ + tol
    regime = 'off';
elseif eps_ > 0 && abs(S) < eps_ - tol
    regime = 'medium';
else
    G = ztl_eom_sun(Y, P, 'off');  sd = dS_dY(Y,P)*G;
    if eps_ == 0
        if real(sd) > 0, regime = 'off'; else, regime = 'on'; end
    else
        if abs(S - eps_) < abs(S + eps_)
            if real(sd) > 0, regime = 'off'; else, regime = 'medium'; end
        else
            if real(sd) < 0, regime = 'on'; else, regime = 'medium'; end
        end
    end
end
end

function newRegime = next_regime(regime, ie, eps_)
switch regime
    case 'on',  newRegime = tern(eps_ > 0, 'medium', 'off');
    case 'off', newRegime = tern(eps_ > 0, 'medium', 'on');
    case 'medium'
        if ie == 1, newRegime = 'off'; else, newRegime = 'on'; end
end
end

function row = dS_dY(Y, P)
% dS/dY [1x15] for S = 1 - ||lam_v|| c/m - lam_m (t-component is 0).
m = Y(7);  lam_v = Y(11:13);  lamvMag = sqrt(sum(lam_v.^2));
row = zeros(1,15);
row(7)     =  P.c*lamvMag/m^2;
row(11:13) = -(P.c/m)*(lam_v.'/lamvMag);
row(14)    = -1;
end

function v = getdef(s,f,d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function v = tern(c,a,b)
if c, v = a; else, v = b; end
end
