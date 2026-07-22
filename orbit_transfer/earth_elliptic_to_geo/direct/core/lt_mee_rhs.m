function [dXdL, Ldot] = lt_mee_rhs(X, U, par)
% LT_MEE_RHS  Gauss variational equations for low-thrust 2-body motion in
% Modified Equinoctial Elements (MEE), independent variable L (true longitude).
%
% State X = [P; ex; ey; hx; hy; m; t] evolves per the paper's Gauss-equation
% block (RTN thrust components), converted from d/dt to d/dL via Ldot. Written
% without norm/abs/max/if on state-dependent quantities and with trig taken on
% an inline floor-based wrap of L into [0,2*pi) -- MATLAB's own documented
% definition mod(a,m) = a - m.*floor(a./m), value-identical to mod(L,2*pi) for
% all real L -- so it evaluates on both numeric doubles and CasADi MX
% symbolics (MX has no overload of the builtin mod(), but floor() works).
%
% INPUTS:
%   X   - State [P; ex; ey; hx; hy; m; t] [7x1]
%   U   - Control [beta(3); thr] with beta a unit RTN thrust direction
%         (radial, transverse, normal) and thr in [0,1] throttle [4x1]
%   par - Struct from kepler_lt_params, with extra field par.L set to the
%         independent-variable value (true longitude) at this node [struct].
%         Optional par.pert (from lunar_params): opt-in lunar third-body
%         term (direct + indirect), added to the RTN force of every Gauss
%         equation including Ldot. Absent, empty, or par.pert.gain == 0
%         takes the untouched nominal branch (bitwise-identical to the
%         pre-2026-07-22 code).
%
% OUTPUTS:
%   dXdL - d/dL of state = (dX/dt)/Ldot [7x1]
%   Ldot - dL/dt [scalar]
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004, p.6 (Gauss equations).
%   [2] earth_elliptic_to_geo/process/DESIGN.md sec 2 (problem statement).
%   [3] docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md
%       (lunar third-body term: constants, RTN frame, gain continuation).

P  = X(1);  ex = X(2);  ey = X(3);  hx = X(4);  hy = X(5);  m = X(6);
L  = par.L;  Tm = par.Tmax;  c = par.c;  mu = par.mu;

thr = U(4);
q = thr*U(1);  s = thr*U(2);  w = thr*U(3);

Lw = L - 2*pi*floor(L/(2*pi));
cL = cos(Lw);  sL = sin(Lw);
Z  = 1 + ex*cL + ey*sL;
A1 = ex + (1+Z)*cL;
A2 = ey + (1+Z)*sL;
Xh = 1 + hx^2 + hy^2;
hterm = hx*sL - hy*cL;
sqPmu = sqrt(P/mu);

% Opt-in lunar third-body term (spec D1/D3; review amendment A 2026-07-22):
% pertOn is a PLAIN-DOUBLE decision at graph-build time, so the nominal
% branch is the literal pre-2026-07-22 code -- bitwise-identical when pert
% is absent OR gain == 0.
pertOn = isfield(par,'pert') && ~isempty(par.pert) && par.pert.gain > 0;
if ~pertOn
    Pdot  = (2*Tm/m)*(P*sqPmu) * (s/Z);
    exdot = (Tm/m)*sqPmu*(1/Z)*( Z*sL*q + A1*s - ey*hterm*w );
    eydot = (Tm/m)*sqPmu*(1/Z)*(-Z*cL*q + A2*s + ex*hterm*w );
    hxdot = (Tm/(2*m))*sqPmu*(Xh/Z)*cL*w;
    hydot = (Tm/(2*m))*sqPmu*(Xh/Z)*sL*w;
    % NOTE: the paper's printed L-dot equation (p.6) omits Tmax on the thrust
    % term -- a typo. The same page's compact form xdot = a(x) + (Tmax/m)*B(x)*u,
    % and the standard Walker/Betts MEE Gauss equations, both carry Tmax/m on
    % every thrust term (compare hxdot/hydot above, which correctly have
    % Tm/(2*m)). Fixed here to (Tm/m).
    Ldot  = sqrt(mu/P^3)*Z^2 + (Tm/m)*sqPmu*(1/Z)*hterm*w;
else
    % Total RTN specific force: thrust + lunar third body (direct + MANDATORY
    % indirect term; pert.gain = the mu-continuation knob; PURE acceleration,
    % no mdot coupling).
    pM = par.pert;
    r  = P/Z;  alpha2 = hx^2 - hy^2;
    rx = (r/Xh)*(cL + alpha2*cL + 2*hx*hy*sL);
    ry = (r/Xh)*(sL - alpha2*sL + 2*hx*hy*cL);
    rz = (2*r/Xh)*(hx*sL - hy*cL);
    Rx = rx/r;  Ry = ry/r;  Rz = rz/r;
    Nx = 2*hy/Xh;  Ny = -2*hx/Xh;  Nz = (1 - hx^2 - hy^2)/Xh;
    Tx = Ny*Rz - Nz*Ry;  Ty = Nz*Rx - Nx*Rz;  Tz = Nx*Ry - Ny*Rx;
    tState = X(7);
    ang = pM.nM*tState + pM.phi0;
    rMx = pM.DM*cos(ang);  rMy = pM.DM*sin(ang);          % Moon in ref plane
    dx = rMx - rx;  dy = rMy - ry;  dz = -rz;
    d3  = (dx^2 + dy^2 + dz^2 + 1e-12)^1.5;               % sep >= 8 LU; guard inert
    DM3 = pM.DM^3;
    gm  = pM.gain * pM.muM;
    aX = gm*(dx/d3 - rMx/DM3);                            % direct + indirect
    aY = gm*(dy/d3 - rMy/DM3);
    aZ = gm*(dz/d3);                                      % Moon z == 0
    fR = (Tm/m)*q + (Rx*aX + Ry*aY + Rz*aZ);
    fT = (Tm/m)*s + (Tx*aX + Ty*aY + Tz*aZ);
    fN = (Tm/m)*w + (Nx*aX + Ny*aY + Nz*aZ);
    Pdot  = 2*(P*sqPmu) * (fT/Z);
    exdot = sqPmu*(1/Z)*( Z*sL*fR + A1*fT - ey*hterm*fN );
    eydot = sqPmu*(1/Z)*(-Z*cL*fR + A2*fT + ex*hterm*fN );
    hxdot = (1/2)*sqPmu*(Xh/Z)*cL*fN;
    hydot = (1/2)*sqPmu*(Xh/Z)*sL*fN;
    Ldot  = sqrt(mu/P^3)*Z^2 + sqPmu*(1/Z)*hterm*fN;
end
mdot  = -(Tm/c)*thr;             % thrust only: gravity costs no propellant
tdot  = 1;

dXdt  = [Pdot; exdot; eydot; hxdot; hydot; mdot; tdot];
% Guard the d/dt -> d/dL division against IPOPT trial steps with Ldot <= 0.
% The NLP's `Ldot >= LdotMin` path constraint only holds at FEASIBLE points, but
% CasADi evaluates this graph during line searches BEFORE the constraint bites,
% so an unguarded 1/Ldot yields NaN -> Invalid_Number_Detected (external review,
% GPT-5.6 + Gemini, 2026-07-19). LdotFloor (default 1e-6) is far below LdotMin
% (default 1e-3), so the guard is INERT at every feasible point and does not
% change certified results -- it only prevents a NaN at pathological trial
% points. (Numeric re-check calls pass doubles; MX build passes symbolics.)
LdotFloor = 1e-6;
if isfield(par, 'LdotFloor') && ~isempty(par.LdotFloor), LdotFloor = par.LdotFloor; end
if isnumeric(Ldot), Ldiv = max(Ldot, LdotFloor); else, Ldiv = fmax(Ldot, LdotFloor); end
dXdL  = dXdt / Ldiv;
end
