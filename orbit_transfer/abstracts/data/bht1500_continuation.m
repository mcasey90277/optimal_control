function out = bht1500_continuation(opts)
% BHT1500_CONTINUATION  Minimum-time DRO->tulip transfer for the Busek BHT-1500.
%
% Reproduces every number in abstract_DRO2tulip_mintime_v3.txt and in
% BHT1500_DRO2tulip_solution.md. Run it and the values should come back
% identical; if they do not, something in pumpkyn/pumpkynPie has moved and the
% abstract needs re-checking before submission.
%
% WHAT IT DOES. Continues the converged minimum-time solution from Darin's
% demo (demos/lowThrustDRO2Tulip.m) to the BHT-1500 operating point of 101 mN
% and 1710 s, walking thrust and specific impulse together in adaptive steps
% and re-seeding each solve from the previous converged costates.
%
% WHY IT STARTS AT 0.1 N, NOT 0.07 N. The demo's ACTIVE case is 70 mN with
% Isp 900 s, and continuing from there nearly stalls: ode45 fails at
% tau ~ 4.06 with the step size driven below 1.4e-14 -- the grazing-singularity
% mode that shows up as revolution count grows -- and the continuation crawled
% to s = 0.083 after repeated step halvings. The demo's THREE COMMENTED-OUT
% costate blocks are in fact converged solutions at its three commented thrust
% tiers: their final times (2.709 / 4.015 / 6.189) line up with 0.1 / 0.07 /
% 0.055 N. Seeding from the 0.1 N block makes the continuation almost pure Isp,
% on a trajectory 32% shorter and far better conditioned. From there it takes
% four steps with zero failures. For this problem the seed's thrust tier
% matters more than the step schedule.
%
% HOW CONVERGENCE IS JUDGED. pumpkyn.cr3bp.tfMin returns no success flag, so
% each step is verified independently by re-propagating and checking the four
% conditions the solver enforces: r(tf) = r_f, v(tf) = v_f, lambda_m(tf) = 0
% and H(tf) = 0. Note this re-integration uses ode45 at default tolerances and
% is STRICTER than the solver's own residual -- the published 70 mN baseline
% re-propagates to 5.0e-06 where fsolve reports ~1e-11. The accepted BHT-1500
% steps all land at ~1.5e-09, tight by either standard.
%
% PREREQUISITE: pumpkynPie and pumpkyn, both at or after the 2026-07-30 pull
% (pumpkynPie 47be599, pumpkyn 5f5ca31). They must be updated TOGETHER -- the
% newer SatelliteAnimator passes 'AddLighting' to pumpkyn's showMoon, so
% pulling one alone breaks the other.
%
% INPUTS:
%   opts - struct (optional):
%          .thrustN   target thrust, N            [default 0.101]
%          .ispS      target specific impulse, s  [default 1710]
%          .m0kg      wet mass, kg                [default 150]
%          .pieRoot   pumpkynPie directory
%                     [default ~/Desktop/proj7/external/pumpkynPie]
%          .tol       accepted terminal residual  [default 1e-8]
%          .verbose                               [default true]
%
% OUTPUTS:
%   out - struct: .lambda0 [7x1] converged costates, .tf_ND, .tf_days,
%         .dV_kms, .propKg, .mfKg, .residual, .steps (the continuation table)
%
% REFERENCES:
%   [1] demos/lowThrustDRO2Tulip.m (the case this continues from).
%   [2] Busek Hall thruster line, https://www.busek.com/hall-thrusters --
%       BHT-1500 published at 101 mN, 1710 s, 1500 W, measured on xenon.
%   [3] Zhang, Topputo, Bernelli-Zazzera & Zhao, "Low-Thrust Minimum-Fuel
%       Optimization in the Circular Restricted Three-Body Problem," JGCD 2015,
%       doi:10.2514/1.G001080.

if nargin < 1, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
T1      = d('thrustN', 0.101);
I1      = d('ispS',    1710);
m0      = d('m0kg',    150);
pieRoot = d('pieRoot', fullfile(getenv('HOME'),'Desktop','proj7','external','pumpkynPie'));
tol     = d('tol',     1e-8);
verbose = d('verbose', true);

here = pwd;  cleaner = onCleanup(@() cd(here)); %#ok<NASGU>
cd(pieRoot);  startup();

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;
g0     = 9.80665*tStar^2/(1000*lStar);

% --- endpoints, built exactly as the demo builds them -----------------------
tauND0 = 1.0;
[~,rvND0] = pumpkynPie.cr3bp.getDRO(tauND0);
rvND0 = pumpkyn.cr3bp.cont_np(rvND0, tauND0, muStar, 1e-12);
[~,rvNDv] = pumpkyn.cr3bp.prop(tauND0, rvND0, muStar);

Np = 7;  tauNDf = 5*2*pi/6;
[~,rvNDF] = pumpkyn.cr3bp.getTulip(tauNDf, Np, -1);
rvNDF = pumpkyn.cr3bp.cont_np(rvNDF, tauNDf, muStar, 1e-12);
[~,rvND] = pumpkyn.cr3bp.prop(tauNDf, rvNDF, muStar);
dvTheta = pumpkyn.util.bsxAng(rvND(:,4:6), rvNDv(1,4:6), 2);
[~,idx] = max(dvTheta);
rvNDF = rvND(idx,:);

ndT = @(TN)   (TN/m0)*tStar^2/(lStar*1000);
ndC = @(IspS) (IspS/tStar)*g0;

% --- the 0.1 N seed: the demo's FIRST commented-out costate block ------------
lam = [-5.92752264703522;  1.65240794854372; 10.9297048382824; ...
       -0.308361759434294; 3.07442715357052; -1.05073539346021; ...
        2.30423587936813;  2.70874702871685];
T0 = 0.100;  I0 = 900;

steps = struct('s',{},'thrustN',{},'ispS',{},'residual',{},'tf_days',{}, ...
               'dV_kms',{},'propKg',{});
if verbose
    fprintf('\n%-6s %8s %8s | %11s %10s %9s %9s\n', ...
            's','T(mN)','Isp(s)','residual','tf(days)','dV(km/s)','prop(kg)');
end

s = 0;  step = 0.25;  nfail = 0;
[res0,~,rv0] = local_check(lam, ndT(T0), ndC(I0), rvND0, rvNDF, muStar);
if verbose
    fprintf('%-6.3f %8.1f %8.0f | %11.2e %10.3f %9.4f %9.3f  SEED\n', 0, T0*1000, I0, ...
        res0, lam(8)*tStar/86400, local_dv(ndC(I0),rv0(end,7),lStar,tStar), m0*(1-rv0(end,7)));
end

while s < 1 - 1e-12
    sn = min(1, s + step);
    Tn = T0 + sn*(T1-T0);
    In = I0 + sn*(I1-I0);
    res = inf;
    try
        solN = pumpkyn.cr3bp.tfMin(rvND0, rvNDF, lam, ndT(Tn), ndC(In), muStar);
        [res,~,rvN] = local_check(solN, ndT(Tn), ndC(In), rvND0, rvNDF, muStar);
    catch
        res = inf;
    end
    if res < tol
        s = sn;  lam = solN;
        dv = local_dv(ndC(In), rvN(end,7), lStar, tStar);
        steps(end+1) = struct('s',s,'thrustN',Tn,'ispS',In,'residual',res, ...
            'tf_days',lam(8)*tStar/86400,'dV_kms',dv,'propKg',m0*(1-rvN(end,7))); %#ok<AGROW>
        if verbose
            fprintf('%-6.3f %8.1f %8.0f | %11.2e %10.3f %9.4f %9.3f\n', s, Tn*1000, In, ...
                    res, steps(end).tf_days, dv, steps(end).propKg);
        end
        step = min(0.25, step*1.5);
    else
        step = step/2;  nfail = nfail + 1;
        if verbose
            fprintf('%-6s %8.1f %8.0f | %11.2e  step halved -> %.4f\n', ...
                    '(fail)', Tn*1000, In, res, step);
        end
        if step < 1e-3
            error('bht1500_continuation:stalled', ...
                  'continuation stalled at s = %.4f (residual %.2e)', s, res);
        end
    end
end

[res,~,rv] = local_check(lam, ndT(T1), ndC(I1), rvND0, rvNDF, muStar);
mf = rv(end,7);
out = struct('lambda0', lam(1:7), 'tf_ND', lam(8), 'tf_days', lam(8)*tStar/86400, ...
    'dV_kms', local_dv(ndC(I1), mf, lStar, tStar), 'propKg', m0*(1-mf), ...
    'mfKg', m0*mf, 'residual', res, 'steps', steps, 'nFail', nfail);

if verbose
    fprintf('\n===== BHT-1500: %g mN, Isp %g s, %g kg =====\n', T1*1000, I1, m0);
    fprintf('  transfer   : %.10f ND = %.3f days\n', out.tf_ND, out.tf_days);
    fprintf('  delta-V    : %.4f km/s\n', out.dV_kms);
    fprintf('  propellant : %.3f kg  (%.2f%% of wet mass)\n', out.propKg, 100*out.propKg/m0);
    fprintf('  final mass : %.3f kg\n', out.mfKg);
    fprintf('  residual   : %.2e   (%d step failures)\n', out.residual, nfail);
end
end

% ---------------------------------------------------------------------------
function [res, tau, rv] = local_check(sol, TmaxND, cND, rvND0, rvNDF, muStar)
% LOCAL_CHECK  Independent verification of one solve.
% Re-propagates and checks the four terminal conditions tfMin enforces, since
% tfMin itself returns no success flag.
% INPUTS:  sol [8x1]; TmaxND; cND; rvND0 [1x6]; rvNDF [1x6]; muStar
% OUTPUTS: res - worst terminal residual; tau, rv - the propagated arc
[tau, rv] = pumpkyn.cr3bp.tfMinProp(sol(8), [rvND0, 1, sol(1:7)'], TmaxND, cND, muStar);
[~, H] = pumpkyn.cr3bp.tfMinEoM(tau(end), rv(end,:)', TmaxND, cND, muStar);
res = max([norm(rv(end,1:6) - rvNDF), abs(rv(end,14)), abs(H)]);
end

% ---------------------------------------------------------------------------
function dv = local_dv(cND, mfFrac, lStar, tStar)
% LOCAL_DV  Dimensional delta-V from the ND exhaust velocity and mass fraction.
% INPUTS: cND; mfFrac - final mass fraction; lStar [km]; tStar [s]
% OUTPUTS: dv [km/s]
dv = (cND*lStar/tStar) * log(1/mfFrac);
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
