function R = run_dro_tulip(opts)
% RUN_DRO_TULIP  Front door: reproduces every number in ../FINDINGS.md.
%
% One call regenerates the whole record -- the reference characterization, the
% mesh sweep, and the lunar-approach table that shows the problem is ill-posed
% without a minimum-altitude constraint.
%
% STAGES
%   1. ENDPOINTS   from pumpkyn's catalogs (getter + cont_np), identical to the
%                  construction in pumpkynPie demos/lowThrustDRO2Tulip.m so the
%                  direct and indirect formulations solve the SAME problem.
%   2. REFERENCE   the indirect min-time solution (pumpkyn.cr3bp.tfMin) from the
%                  demo's own converged costates, characterized: revolutions,
%                  Earth and Moon distances, mass.
%   3. SWEEP       the direct min-time twin at several mesh densities, each
%                  reporting t_f, feasibility, throttle saturation and closest
%                  lunar approach.
%
% WHAT THE OUTPUT SHOWS. Solutions are machine-tight yet t_f depends
% non-monotonically on N, and the fastest undercuts the indirect reference by
% flying a much deeper lunar flyby. See ../FINDINGS.md.
%
% PREREQUISITES: pumpkynPie + pumpkyn (pull them TOGETHER), CasADi 3.7.0.
%
% INPUTS:
%   opts - struct (optional):
%          .N        mesh densities to sweep       [default [400 800 1600]]
%          .thrustN  thrust, N                     [default 0.07]
%          .ispS     specific impulse, s           [default 900]
%          .m0kg     wet mass, kg                  [default 150]
%          .save     write results/*.mat           [default true]
%          .verbose                                [default true]
%
% OUTPUTS:
%   R - struct: .rv0 .rvf .p .ref (reference characterization)
%       .sweep (struct array per N: .N .tf .maxDefect .thrMin .minR2 .minR1
%       .altKm .status) .Tmax .c
%
% REFERENCES:
%   [1] ../FINDINGS.md — the record this reproduces.
%   [2] pumpkynPie demos/lowThrustDRO2Tulip.m — the indirect reference.

if nargin < 1, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
Ns      = d('N', [400 800 1600]);
TmaxN   = d('thrustN', 0.07);
IspS    = d('ispS', 900);
m0      = d('m0kg', 150);
doSave  = d('save', true);
verbose = d('verbose', true);

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'lib'));
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
cwd0 = pwd;  cleaner = onCleanup(@() cd(cwd0)); %#ok<NASGU>
cd(fullfile(getenv('HOME'),'Desktop','proj7','external','pumpkynPie'));  startup();
cd(cwd0);

%% 1. endpoints
[rv0, rvf, p] = dro_tulip_endpoints();
g0   = 9.80665*p.tStar^2/(1000*p.lStar);
Tmax = (TmaxN/m0)*p.tStar^2/(p.lStar*1000);
c    = (IspS/p.tStar)*g0;
mu   = p.muStar;

%% 2. the indirect reference, from the demo's converged costates
lam = [13.9579470518969; 6.26766054139842; 7.02087277890196; -0.969930036390389; ...
       -5.04449989501292; -3.11628049459521; 5.50150618633038; 4.01524259262941];
sol = pumpkyn.cr3bp.tfMin(rv0, rvf, lam, Tmax, c, mu);
[tauR, rvR] = pumpkyn.cr3bp.tfMinProp(sol(8), [rv0, 1, sol(1:7)'], Tmax, c, mu);
r1R = vecnorm(rvR(:,1:3) - [-mu 0 0], 2, 2);
r2R = vecnorm(rvR(:,1:3) - [1-mu 0 0], 2, 2);
thR = unwrap(atan2(rvR(:,2), rvR(:,1)-(1-mu)));
ref = struct('tf', sol(8), 'tf_days', sol(8)*p.tStar/86400, ...
    'revsMoon', (max(thR)-min(thR))/(2*pi), 'minR1', min(r1R), 'maxR1', max(r1R), ...
    'minR2', min(r2R), 'altKm', min(r2R)*p.lStar - 1737, 'mf', rvR(end,7), ...
    'costates', sol(1:7));

if verbose
    fprintf('\n===== REFERENCE (indirect, pumpkyn.cr3bp.tfMin) =====\n');
    fprintf('  t_f            : %.6f ND = %.3f days\n', ref.tf, ref.tf_days);
    fprintf('  revs about Moon: %.2f      <- why no Sundman is needed\n', ref.revsMoon);
    fprintf('  Earth distance : %.4f .. %.4f (never approaches Earth)\n', ref.minR1, ref.maxR1);
    fprintf('  closest Moon   : %.5f ND = %.0f km altitude\n', ref.minR2, ref.altKm);
    fprintf('  mass fraction  : 1.0 -> %.6f\n', ref.mf);
end

%% 3. the direct twin, swept over mesh density
sweep = struct('N',{},'tf',{},'maxDefect',{},'thrMin',{},'minR2',{},'minR1',{}, ...
               'altKm',{},'status',{},'wall',{});
if verbose
    fprintf('\n===== DIRECT TWIN (casadi_mintime_dro), crude seed =====\n');
    fprintf('%-7s %10s %10s %10s %12s %13s %s\n', ...
        'N','t_f','vs ref','defect','min thr','Moon alt km','status');
end
for k = 1:numel(Ns)
    N = Ns(k);  t0 = tic;
    o = casadi_mintime_dro(rv0, rvf, Tmax, c, mu, N, [], [], [], struct('maxIter',3000));
    w = toc(t0);
    r2 = vecnorm(o.X(1:3,:) - [1-mu;0;0], 2, 1);
    r1 = vecnorm(o.X(1:3,:) - [-mu;0;0], 2, 1);
    sweep(k) = struct('N',N,'tf',o.tf,'maxDefect',o.maxDefect,'thrMin',o.thrMin, ...
        'minR2',min(r2),'minR1',min(r1),'altKm',min(r2)*p.lStar-1737, ...
        'status',o.ipoptStatus,'wall',w); %#ok<AGROW>
    if verbose
        fprintf('%-7d %10.6f %+10.4f %10.1e %10.6f %13.0f %s\n', ...
            N, o.tf, o.tf-ref.tf, o.maxDefect, o.thrMin, sweep(k).altKm, o.ipoptStatus);
    end
    if doSave
        rd = fullfile(here,'results');  if ~isfolder(rd), mkdir(rd); end
        save(fullfile(rd, sprintf('mintime_N%d.mat', N)), 'o','rv0','rvf','p','Tmax','c','-v7.3');
    end
end

if verbose
    [~,ix] = min([sweep.tf]);
    fprintf('\n  FASTEST: N = %d at t_f = %.6f, %.1f%% below the reference,\n', ...
        sweep(ix).N, sweep(ix).tf, 100*(ref.tf-sweep(ix).tf)/ref.tf);
    fprintf('  flying a %.0f km lunar pass against the reference''s %.0f km.\n', ...
        sweep(ix).altKm, ref.altKm);
    fprintf('  A minimum-time problem cannot be beaten -- so the reference is a\n');
    fprintf('  LOCAL minimum, and without an altitude constraint the problem is\n');
    fprintf('  ill-posed. See ../FINDINGS.md.\n');
end

R = struct('rv0',rv0,'rvf',rvf,'p',p,'ref',ref,'sweep',sweep,'Tmax',Tmax,'c',c);
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
