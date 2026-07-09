function out = run_sundman_from_seed(N, pSund, seedFile, maxIter)
% RUN_SUNDMAN_FROM_SEED  Sundman-regularized min-fuel solve seeded from a
% certified, collocation-feasible time-mesh solution (not an ODE burn+coast).
%
% The plain burn+coast warm start is only ODE-consistent: after resampling
% onto a uniform-sigma Sundman mesh its terminal node drifts off the target
% over ~40 revolutions, so IPOPT declares local infeasibility. Seeding from a
% solution that ALREADY hits rvf exactly (endpoint err 0, defect 2e-15)
% removes that drift -- the only remaining job for IPOPT is to close the
% modest resampling defect and (later, under tf-continuation) grow switches.
%
% INPUTS:
%   N        - number of Sundman segments [default 2000]
%   pSund    - Sundman power kappa=r1^pSund [default 1.5]
%   seedFile - .mat with .nlp.X (7xM: r,v,m), .nlp.U (4xM: [w;s] cone),
%              .sigma (Mx1), .tf, .rv0, .rvf  [default the energy-seed soln]
%   maxIter  - IPOPT max iterations [default 3000]
%
% OUTPUTS:
%   out - struct from casadi_minfuel_sundman (+ prints dV, switches)

if nargin<1||isempty(N), N=2000; end
if nargin<2||isempty(pSund), pSund=1.5; end
here=fileparts(mfilename('fullpath'));
if nargin<3||isempty(seedFile), seedFile=fullfile(here,'minfuel_from_energy_seed.mat'); end
if nargin<4||isempty(maxIter), maxIter=3000; end
addpath(here); run(fullfile(here,'setup_paths.m'));

muStar=0.012150585609624; lStar=389703.264829278; tStar=382981.289129055;
m0kg=15; g0=9.80665*tStar^2/(1000*lStar);
Tmax=(0.025/m0kg)*tStar^2/(lStar*1000); c=(2100/tStar)*g0;

S = load(seedFile);
Xs = S.nlp.X;                 % 7 x M  (r;v;m)
Us = S.nlp.U;                 % 4 x M  ([w;s] cone form)
tf = S.tf;  rv0 = S.rv0;  rvf = S.rvf;
sg = S.sigma(:);  sg = (sg-sg(1))/(sg(end)-sg(1));   % normalize 0..1
tSeed = sg*tf;                % uniform-time mesh -> physical time per node

% cone form [w;s] -> unit direction alpha + throttle s
wcol = Us(1:3,:);  s_seed = Us(4,:);
alpha = wcol ./ max(sqrt(sum(wcol.^2,1)),1e-9);

% --- map time -> Sundman tau: dtau = dt/kappa, kappa=r1^pSund ---
r1 = sqrt((Xs(1,:)+muStar).^2 + Xs(2,:).^2 + Xs(3,:).^2).';
kap = r1.^pSund;
dt = diff(tSeed);
dtau = dt.*0.5.*(1./kap(1:end-1)+1./kap(2:end));
tau = [0; cumsum(dtau)]; tauf0 = tau(end); sig_i = tau/tauf0;

% NO-RESAMPLE: use the seed's OWN nodes (sigma = tau/tauf0). Downsampling a
% 40-rev oscillatory trajectory onto a uniform-sigma mesh leaves an
% irreducible ~1e-2 defect that pins IPOPT in restoration; the seed's own
% nodes make the only initial infeasibility the small time-trap vs
% Sundman-trap mismatch, which IPOPT closes to ~1e-14. N is ignored (kept for
% signature compatibility). This matches run_sundman_homotopy.
[sigma,ku] = unique(sig_i,'stable');
X0 = [Xs(:,ku); tSeed(ku).'];
U0 = [alpha(:,ku); s_seed(ku)];
N  = numel(sigma)-1; %#ok<NASGU>

% pin endpoints exactly
X0(1:6,1)=rv0(:); X0(7,1)=1; X0(8,1)=0; X0(1:6,end)=rvf(:); X0(8,end)=tf;

fprintf('SUNDMAN-from-seed: N=%d pSund=%.2f tf=%.4f tauf0=%.4g (seed switches=%d)\n', ...
        N, pSund, tf, tauf0, sum(abs(diff(s_seed>0.5))));
out = casadi_minfuel_sundman(sigma, tf, rv0, rvf, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter);
dV = c*log(1/out.mf)*lStar/tStar;
fprintf('\n=== SUNDMAN-from-seed RESULT (pSund=%.2f) ===\n', pSund);
fprintf('success=%d status=%s\n', out.success, out.ipoptStatus);
fprintf('maxDefect=%.2g maxUnit=%.2g tauf=%.4g\n', out.maxDefect, out.maxUnit, out.tauf);
fprintf('switches=%d bang-bang=%.1f%% prop=%.4f kg dV=%.4f km/s\n', ...
        out.switches, 100*out.edge, m0kg*(1-out.mf), dV);
fprintf('SUNDMAN_PASS=%d\n', out.maxDefect<1e-6 && out.success);
end
