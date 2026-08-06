function B = survey_period_bounds(outMat)
%% Purpose:
%
%   Finds the "REASONABLE" period ranges for the DRO and tulip families, per
%   Darin's criteria: the orbit must not approach the Moon dangerously, and
%   should stay within 100 Mm (100,000 km) of the Moon -- near enough to be
%   lunar-vicinity, far enough out to take in L1 and L2.
%
%   Pure propagation, no optimization: sweep each family's parameter, build
%   the orbit with pumpkyn, and measure
%
%     periAltKm   closest approach to the lunar SURFACE (altitude, km)
%     maxMoonKm   farthest distance from the Moon (km)
%     periodDays  the orbit period
%
%   An orbit is ADMISSIBLE when periAltKm >= altFloorKm (default 500, the
%   same floor the transfer campaigns use) and maxMoonKm <= 100,000. For the
%   DRO, tau IS the period in ND time (getDRO selects by period). For the
%   tulip, the family is resonant -- Np petals lock the period to
%   (Np-2)/(Np-1) lunar months -- so the sweep is over petal count Np and
%   branch pm, with tau at each Np's resonant value.
%
%% Inputs:
%
%  outMat                   char                    Output .mat for the
%                                                   survey tables (optional;
%                                                   '' = no file)
%
%% Outputs:
%
%  B                        struct
%   .dro                    table-like struct        tau, periodDays,
%                                                    periAltKm, maxMoonKm,
%                                                    admissible per member
%   .tulip                  table-like struct        Np, pm, tau, periodDays,
%                                                    periAltKm, maxMoonKm,
%                                                    admissible per member
%   .droTauRange            [1 x 2]                  admissible tau (period)
%                                                    bounds for the DRO
%   .tulipNpAdmissible      vector                   petal counts admissible
%                                                    on either branch
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, outMat = ''; end

   muStar = 0.012150585609624;
    lStar = 389703.264829278;            % km
    tStar = 382981.289129055;            % s
  rMoonKm = 1737.4;
altFloorKm = 500;                        % same floor as the transfer work
 maxDistKm = 100000;                     % Darin's 100 Mm lunar-vicinity box
    moonND = [1-muStar, 0, 0];

%% ---- DRO family: tau IS the period (ND). Sweep it. ---------------------
  tauGrid = [0.05:0.05:0.5, 0.6:0.1:2.0, 2.25:0.25:5.0];
      dro = struct('tau',[], 'periodDays',[], 'periAltKm',[], ...
                   'maxMoonKm',[], 'admissible',[]);
fprintf('DRO family (tau = period in ND; 1 ND = %.3f days):\n', tStar/86400);
fprintf('   tau   period(d)   periAlt(km)   maxMoon(km)   admissible\n');
for tau = tauGrid
    try
        [~, rv0] = pumpkynPie.cr3bp.getDRO(tau);
        rv0 = pumpkyn.cr3bp.cont_np(rv0, tau, muStar, 1e-12);
        [tt, rv] = pumpkyn.cr3bp.prop(tau, rv0, muStar);
        dM = vecnorm(rv(:,1:3) - moonND, 2, 2) * lStar;
        pa = min(dM) - rMoonKm;
        mx = max(dM);
        ok = pa >= altFloorKm && mx <= maxDistKm;
        dro.tau(end+1) = tau;                      %#ok<*AGROW>
        dro.periodDays(end+1) = tt(end)*tStar/86400;
        dro.periAltKm(end+1) = pa;
        dro.maxMoonKm(end+1) = mx;
        dro.admissible(end+1) = ok;
        fprintf('  %5.2f   %7.3f   %11.0f   %11.0f   %s\n', ...
                tau, tt(end)*tStar/86400, pa, mx, admstr(ok));
    catch
        fprintf('  %5.2f   (family generator failed)\n', tau);
    end
end

%% ---- Tulip family: resonance locks the period; sweep petals Np, pm -----
% Np petals <-> tau = (Np-2)/(Np-1) * 2*pi  (verified: Np=7 -> 5*2*pi/6)
    tulip = struct('Np',[], 'pm',[], 'tau',[], 'periodDays',[], ...
                   'periAltKm',[], 'maxMoonKm',[], 'admissible',[]);
fprintf('\nTulip family (period locked by resonance to (Np-2)/(Np-1) months):\n');
fprintf('   Np  pm    tau     period(d)   periAlt(km)   maxMoon(km)   admissible\n');
for Np = 3:14
    for pm = [-1 1]
        tau = (Np-2)/(Np-1) * 2*pi;
        try
            [~, rv0] = pumpkyn.cr3bp.getTulip(tau, Np, pm);
            rv0 = pumpkyn.cr3bp.cont_np(rv0, tau, muStar, 1e-12);
            [tt, rv] = pumpkyn.cr3bp.prop(tau, rv0, muStar);
            dM = vecnorm(rv(:,1:3) - moonND, 2, 2) * lStar;
            pa = min(dM) - rMoonKm;
            mx = max(dM);
            ok = pa >= altFloorKm && mx <= maxDistKm;
            tulip.Np(end+1) = Np;  tulip.pm(end+1) = pm;
            tulip.tau(end+1) = tau;
            tulip.periodDays(end+1) = tt(end)*tStar/86400;
            tulip.periAltKm(end+1) = pa;
            tulip.maxMoonKm(end+1) = mx;
            tulip.admissible(end+1) = ok;
            fprintf('  %3d  %2d  %6.3f   %8.3f   %11.0f   %11.0f   %s\n', ...
                    Np, pm, tau, tt(end)*tStar/86400, pa, mx, admstr(ok));
        catch
            fprintf('  %3d  %2d  %6.3f   (family generator failed)\n', Np, pm, tau);
        end
    end
end

%% ---- The admissible box -----------------------------------------------
adm = dro.tau(logical(dro.admissible));
if isempty(adm), B_dro = [NaN NaN]; else, B_dro = [min(adm) max(adm)]; end
NpAdm = unique(tulip.Np(logical(tulip.admissible)));
B = struct('dro',dro, 'tulip',tulip, 'droTauRange',B_dro, ...
           'tulipNpAdmissible',NpAdm, ...
           'criteria', struct('altFloorKm',altFloorKm, 'maxDistKm',maxDistKm));
fprintf('\nADMISSIBLE BOX: DRO tau (= period, ND) in [%.2f, %.2f];  tulip Np in {%s}\n', ...
        B_dro(1), B_dro(2), num2str(NpAdm));
if ~isempty(outMat), save(outMat, '-struct', 'B'); end
end

% ------------------------------------------------------------------------
function s = admstr(t)
% ADMSTR  yes/NO marker.  INPUTS: t logical.  OUTPUTS: s [char].
if t, s = 'yes'; else, s = 'NO'; end
end
