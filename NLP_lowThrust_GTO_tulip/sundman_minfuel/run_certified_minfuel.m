function best = run_certified_minfuel(maxIter)
% RUN_CERTIFIED_MINFUEL  End-to-end reproduction of the certified sharp
% bang-bang minimum-fuel GTO -> tulip transfer, using the modular library.
%
% Pipeline (each step is a reusable library function):
%   1. cr3bp_lt_params      - CR3BP + low-thrust ND constants
%   2. (seed)               - certified energy-seeded time-mesh solution
%   3. sundman_seed_map     - no-resample map into Sundman coordinates
%   4. sundman_homotopy     - guarded energy->fuel sweep (eps 1 -> 0)
%
% Reproduces: 25-switch bang-bang, defect ~1e-14, terminal error 0,
% propellant ~2.2640 kg, dV ~3.3696 km/s. Saves sundman_minfuel_certified.mat.
%
% INPUTS:
%   maxIter - IPOPT max iters per homotopy step [scalar, default 1500]
%
% OUTPUTS:
%   best - certified solver struct (best point of the homotopy)

if nargin < 1 || isempty(maxIter), maxIter = 1500; end
here = fileparts(mfilename('fullpath'));  addpath(here);
pSund = 1.5;

p = cr3bp_lt_params(0.025, 15, 2100);           % 25 mN, 15 kg, Isp 2100 s

% collocation-feasible time-mesh seed (energy-seeded 3-switch min-fuel)
S   = load(fullfile(here, 'minfuel_from_energy_seed.mat'));
rv0 = S.rv0;  rvf = S.rvf;

% optional consistency check against a freshly built endpoint (needs pumpkyn)
try
    run(fullfile(here, 'setup_paths.m'));
    [rv0c, rvfc] = gto_tulip_endpoints(p);
    fprintf('endpoint check vs pumpkyn: |drv0|=%.1e  |drvf|=%.1e\n', ...
            norm(rv0c-rv0), norm(rvfc-rvf));
catch
    fprintf('(pumpkyn unavailable; using seed-stored endpoints)\n');
end

% map seed into Sundman coordinates (no resample) and sweep eps 1 -> 0
[sigma, X0, U0, tauf0] = sundman_seed_map(S.nlp.X, S.nlp.U, S.tf, S.sigma, ...
                                          pSund, p.muStar, rv0, rvf);
epsSched = [1 0.6 0.35 0.2 0.12 0.07 0.04 0.03 0.022 0.016 0.012 0.009 ...
            0.006 0.004 0.0025 0.0015 0.001 0.0005 0];
saveFile = fullfile(here, 'sundman_minfuel_certified.mat');
best = sundman_homotopy(p, rv0, rvf, sigma, X0, U0, tauf0, pSund, ...
                        epsSched, maxIter, saveFile);
end
