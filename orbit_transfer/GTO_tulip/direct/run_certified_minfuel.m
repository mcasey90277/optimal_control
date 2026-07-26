function best = run_certified_minfuel(maxIter, saveFile)
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
%   maxIter  - IPOPT max iters per homotopy step [scalar, default 1500]
%   saveFile - where to write the result [char]. DEFAULT IS THE PUBLISHED
%              REFERENCE lib/sundman_minfuel_certified.mat, which this
%              OVERWRITES. Pass an explicit path when re-solving to VERIFY
%              rather than to republish -- a re-solve is not bit-identical to
%              the published artifact (see the note below), so overwriting
%              silently replaces the campaign's reference with a different
%              extremal.
%
% OUTPUTS:
%   best - certified solver struct (best point of the homotopy)
%
% A RE-SOLVE DOES NOT REPRODUCE THE PUBLISHED ARTIFACT EXACTLY, and that is a
% known property of this problem, not a defect. Measured re-solves:
%   published artifact   25 switches, m_f 0.849066, dV 3.3696 km/s
%   re-solve 2026-07-21  24 switches,               dV 3.3660 km/s
%   re-solve 2026-07-26  24 switches, m_f 0.849279, dV 3.3644 km/s
% Mass agrees to ~0.1% every time; the switch INTEGER is basin-sensitive even
% at fixed mesh. Treat mass/dV as the reproducible quantities and the switch
% count as a band. Open item: TODO C3 (mesh-band study).

if nargin < 1 || isempty(maxIter), maxIter = 1500; end
here = fileparts(mfilename('fullpath'));  addpath(here);
if nargin < 2 || isempty(saveFile)
    saveFile = fullfile(here, 'lib', 'sundman_minfuel_certified.mat');
    warning('run_certified_minfuel:overwritesReference', ...
        ['about to OVERWRITE the published reference %s. A re-solve lands on a ' ...
         'nearby but different extremal (24 vs 25 switches historically), so this ' ...
         'replaces the campaign reference. To verify instead, pass an explicit ' ...
         'saveFile: run_certified_minfuel(1500, ''/tmp/check.mat'')'], saveFile);
end
pSund = 1.5;

p = cr3bp_lt_params(0.025, 15, 2100);           % 25 mN, 15 kg, Isp 2100 s

% collocation-feasible time-mesh seed (energy-seeded 3-switch min-fuel)
S   = load(fullfile(here, 'lib', 'minfuel_from_energy_seed.mat'));
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
best = sundman_homotopy(p, rv0, rvf, sigma, X0, U0, tauf0, pSund, ...
                        epsSched, maxIter, saveFile);
end
