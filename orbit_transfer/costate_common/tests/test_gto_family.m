function ok = test_gto_family()
% TEST_GTO_FAMILY  Unit test for the 'gto' pseudo-family (Stage B spec 3.1).
%
% Five checks: flagship equality (incl. tau/periodND well-formedness and
% strict monotonicity), rotation equivariance, apsidal geometry, and an
% INDEPENDENT non-apsidal (mean anomaly pi/2) cross-check that guards
% against a mean-vs-true anomaly confusion in the case implementation
% (apsidal points cannot detect this bug because E == nu == M there).
%
%% Inputs:
%
%  (none)
%
%% Outputs:
%
%  ok                       logical                 true iff every check
%                                                   passed
%
%% Revision History:
%  M. Casey                                                   (c) 08/26/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                          % costate_common
addpath(fullfile(fileparts(fileparts(fileparts(here))), ...  % repo root
        'orbit_transfer', 'GTO_tulip', 'indirect', 'min_time'));
ok = true;

% 1. FLAGSHIP EQUALITY: orientDeg=-25 at perigee reproduces mintime_params rv0
[rv0Ref, ~, P] = mintime_params();
[tau, rv, info] = get_family_orbit('gto', struct('orientDeg', -25));
ok = check('perigee state == flagship rv0 (1e-10)', ...
           max(abs(rv(1,:) - rv0Ref(:).')) < 1e-10) && ok;
ok = check('tau starts at 0', tau(1) == 0) && ok;
ok = check('periodND positive and ~ Kepler', ...
           info.periodND > 0 && abs(tau(end) - info.periodND) < 1e-9) && ok;
ok = check('tau strictly increasing (locus parameter, NEVER an epoch offset)', ...
           all(diff(tau) > 0)) && ok;

% 2. ROTATION EQUIVARIANCE: rotating the orientation rotates the whole locus
%    (positions about the barycenter-frame z axis THROUGH EARTH at (-mu,0,0):
%    the ellipse is Earth-centered, so rotate the Earth-relative vector)
th = 37 * pi/180;
Rz = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
[~, rvA] = get_family_orbit('gto', struct('orientDeg', 10));
[~, rvB] = get_family_orbit('gto', struct('orientDeg', 10 + 37));
rE = [-P.muStar 0 0];
posErr = 0; velErr = 0;
for k = [1, 7, 13]                                  % spot rows
    pA = (Rz*(rvA(k,1:3) - rE).').' + rE;
    vA = (Rz*rvA(k,4:6).').';
    posErr = max(posErr, max(abs(pA - rvB(k,1:3))));
    velErr = max(velErr, max(abs(vA - rvB(k,4:6))));
end
ok = check(sprintf('rotation equivariance pos %.1e vel %.1e < 1e-10', ...
           posErr, velErr), posErr < 1e-10 && velErr < 1e-10) && ok;

% 3. GEOMETRY: radii from Earth at frac 0 / 0.5 are perigee / apogee
rp_km = (6378+350); ra_km = (6378+35786);
rPer = norm(rv(1,1:3) - rE) * P.lStar;
kAp = round((numel(tau)+1)/2);                      % M = pi at half period
rAp  = norm(rv(kAp,1:3) - rE) * P.lStar;
ok = check(sprintf('perigee radius %.1f ~ %d km', rPer, rp_km), ...
           abs(rPer - rp_km) < 1.0) && ok;
ok = check(sprintf('apogee radius %.1f ~ %d km', rAp, ra_km), ...
           abs(rAp - ra_km) < 1.0) && ok;

% 4. NON-APSIDAL INDEPENDENT CHECK (Gemini, critical): mean anomaly M=pi/2.
%    Apsidal points cannot distinguish a mean-vs-true anomaly bug (E==nu==M
%    there), so this is the actual arbiter. The expected state is built
%    from scratch here -- own Newton solve for E, own perifocal-ellipse
%    formulas, own rotation, own rotating-frame velocity (subtract omega x
%    r, valid because fromPCI is called at t=0 where the CR3BP rotating
%    frame is instantaneously aligned with the inertial frame, per
%    pumpkyn.cr3bp.primary2PosVel: tau=0 => P2 lies on +x) -- none of it
%    calls pumpkyn.cr3bp.orb2eci/fromPCI. It is then compared against the
%    ACTUAL get_family_orbit output (spline-interpolated onto the exact
%    tau for M=pi/2; verified empirically to agree with the closed form to
%    ~5e-12 on this dense M=2001 grid, far inside the 1e-9 gate).
smaKm  = (6378+350 + 6378+35786)/2;                 % default sma_km
eccGto = (35786-350)/(2*smaKm);                     % default ecc
orient = -25;                                       % flagship orientDeg
muE = 6.67384e-20*(1 - P.muStar)*(5.9736E24 + 7.35E22);

Mtarget = pi/2;
E = Mtarget;                                        % Newton: E - e*sinE = M
for it = 1:50
    dE = -(E - eccGto*sin(E) - Mtarget) / (1 - eccGto*cos(E));
    E  = E + dE;
    if abs(dE) < 1e-14, break; end
end
nu = 2*atan2(sqrt(1+eccGto)*sin(E/2), sqrt(1-eccGto)*cos(E/2));

p_slr = smaKm*(1-eccGto^2);
rm    = p_slr/(1+eccGto*cos(nu));
r_pf  = rm*[cos(nu); sin(nu); 0];
v_pf  = sqrt(muE/p_slr)*[-sin(nu); eccGto+cos(nu); 0];
c = cos(orient*pi/180); s = sin(orient*pi/180);
r_in = [c -s 0; s c 0; 0 0 1]*r_pf;
v_in = [c -s 0; s c 0; 0 0 1]*v_pf;
r_nd = r_in/P.lStar;
v_nd = v_in*(P.tStar/P.lStar);
rvIndep = zeros(1,6);
rvIndep(1:3) = (r_nd + [-P.muStar; 0; 0]).';
rvIndep(4:6) = (v_nd - cross([0;0;1], r_nd)).';

tau_target = info.periodND * (Mtarget/(2*pi));
rvInterp = zeros(1,6);
for c6 = 1:6
    rvInterp(c6) = interp1(tau, rv(:,c6), tau_target, 'spline');
end
errNonApsidal = max(abs(rvInterp - rvIndep));
ok = check(sprintf('non-apsidal (M=pi/2) independent agreement %.1e < 1e-9', ...
           errNonApsidal), errNonApsidal < 1e-9) && ok;

fprintf('test_gto_family: %s\n', string(ok));
end

% ------------------------------------------------------------------------
function ok = check(name, cond)
% CHECK  One gate line.  INPUTS: name;cond  OUTPUTS: ok
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end
