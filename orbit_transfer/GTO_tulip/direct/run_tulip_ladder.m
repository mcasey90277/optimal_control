% RUN_TULIP_LADDER  Front door: walk the GTO->tulip THRUST ladder.
%
% Set the parameters in section 1, run, and get a chained sequence of certified
% min-fuel solutions from thrustStart down to thrustStop, each gated and saved.
%
%   1. PARAMETERS   - thrust range and step, t_f factor, epsilon endpoint
%   2. CHAIN        - eps=1 energy continuation, small steps, rung -> rung
%   3. SHARPEN      - eps: 1 -> epsMin on the requested rung(s)
%   4. GATE         - defect, bound saturation, fingerprint, primer alignment
%   5. REPORT       - one table, plus a per-rung .mat under results/ladder/
%
% WHY A SEPARATE FRONT DOOR FROM run_gto_tulip. run_gto_tulip solves ONE
% transfer at the campaign's nominal 25 mN and deliberately does not expose
% thrust. This script is the off-nominal path, and it exists because the ladder
% turned out to be reachable after all (2026-07-27): 25 -> 20 mN in ~4% steps,
% every rung machine-tight, on an item recorded as blocked by a "topology wall"
% since 2026-07-21. That wall was an artifact of taking one 20% continuation
% step instead of five small ones.
%
% STEP SIZE IS LOAD-BEARING, NOT COSMETIC. A single 25 -> 20 mN jump fails
% outright (inf_du diverges past 1e14, then a MEX crash). ~4% steps converge to
% ~1e-14 at every rung. If you widen stepFrac and rungs start failing, suspect
% the step before suspecting the physics.
%
% WHY THIS USES THE FREE-TIME SOLVER. Below ~21 mN the campaign's own
% casadi_minfuel_sundman (fixed tau_f) stops converging; measured 2026-07-27, it
% walks to 22 mN and fails at 21. casadi_energy_freetf carries Betts' cScale
% slack state (dt/dtau = c*kappa, dc/dtau = 0), which reaches 20 mN. Run with
% moonZone <= 0 it is otherwise IDENTICAL to the tulip formulation -- same
% single-primary Sundman clock, same dynamics, same endpoints. It lives in
% GTO_ELFO because that campaign needed free time first; the dependency is
% declared explicitly below rather than hidden in setup_paths.
%
% KNOWN CEILING (measured 2026-07-27, three hypotheses tested and eliminated):
%   ~20 mN for a min-fuel (eps=0) solution; ~19.5 mN for energy (eps=1) only.
%   NOT resolution   -- 190 nodes/rev at the ceiling, generous.
%   NOT step size    -- 2.5% steps stop in the same place.
%   NOT the t_f rule -- +9% transfer time changes nothing.
%   Both solvers wall within 8% of each other (21 mN fixed, 19.5 mN free),
%   which points at the problem near 19-21 mN rather than the transcription.
%
% OUTPUTS (under results/ladder/):
%   tulip_T<tag>_f<factor>_eps<tag>.mat  - per rung: out, sigma, tauf0, rv0,
%                                          rvf, T, fp, gate
%   ladder_summary.txt                   - the table this prints
%
% REFERENCES:
%   [1] ../TODO.md "Free-span reformulation" -- the full measurement record.
%   [2] doc/gto_tulip_guide.pdf sec 3.4 -- why the mesh is inherited, not chosen.
%   [3] run_gto_tulip.m -- the single-transfer front door at nominal thrust.

%% ------------------------------------------------------------------------
%% 1. PARAMETERS  (edit this section only)
%% ------------------------------------------------------------------------
thrustStart = 0.025;     % [N] start of the chain (25 mN = the certified rung)
thrustStop  = 0.020;     % [N] stop here. Below ~20 mN does not converge -- see
                         %     the KNOWN CEILING note above before lowering it
stepFrac    = 0.04;      % fractional thrust step per rung (~4%; see the note --
                         %     20% fails outright, 2.5% buys one extra rung)
factor      = 1.15;      % t_f / t_f,min, held across rungs (t_f,min ~ 1/T)
epsMin      = 0;         % homotopy endpoint: 1 = energy, 0 = min-fuel
sharpenAll  = false;     % false: sharpen only the FINAL rung (fast)
                         % true : sharpen every rung (a full certified ladder)
doGate      = true;      % run the gate battery on each sharpened rung
rerun       = false;     % true -> ignore existing rung artifacts and redo

% Programmatic override hook (batch/shell drivers), house pattern.
if exist('LADDER_OVERRIDES','var') && isstruct(LADDER_OVERRIDES)
    lovf = fieldnames(LADDER_OVERRIDES);
    for lovk = 1:numel(lovf), eval([lovf{lovk} ' = LADDER_OVERRIDES.(lovf{lovk});']); end
    clear lovf lovk
end

%% ------------------------------------------------------------------------
%% 2. SETUP
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths;

% The free-time solver. Declared here, deliberately, rather than buried in
% setup_paths: it is a real cross-campaign dependency and should be visible to
% anyone reading this file. See "WHY THIS USES THE FREE-TIME SOLVER" above.
elfoDir = fullfile(here, '..', '..', 'GTO_ELFO', 'direct', 'elfo');
assert(isfolder(elfoDir), 'run_tulip_ladder:noELFO', ...
    'free-time solver not found at %s', elfoDir);
addpath(elfoDir);

cfgNom = minfuel_config();
[rv0, rvf] = insertion_states('tulip', 'campaign');
outDir = fullfile(here, 'results', 'ladder');
if ~exist(outDir, 'dir'), mkdir(outDir); end

assert(thrustStop <= thrustStart, 'run_tulip_ladder:range', ...
    'thrustStop (%.4g) must be <= thrustStart (%.4g)', thrustStop, thrustStart);
assert(stepFrac > 0 && stepFrac < 0.25, 'run_tulip_ladder:step', ...
    ['stepFrac = %.3g is outside (0, 0.25). A 20%% step is MEASURED to fail ' ...
     'outright on this ladder; ~4%% is the working value.'], stepFrac);

% The rung sequence: geometric in thrust, always ending exactly on thrustStop.
rungs = thrustStart;
while rungs(end) > thrustStop * (1 + 1e-9)
    rungs(end+1) = max(rungs(end) * (1 - stepFrac), thrustStop);
end

tfmin = @(T) cfgNom.tfMin * (cfgNom.thrustN / T);   % t_f,min ~ 1/T
fprintf(['\n=== TULIP THRUST LADDER ==============================\n' ...
         '  %.1f -> %.1f mN in %d rungs (%.1f%% steps)\n' ...
         '  factor %.3f, epsMin %g, sharpen %s\n' ...
         '======================================================\n'], ...
        thrustStart*1e3, thrustStop*1e3, numel(rungs)-1, 100*stepFrac, ...
        factor, epsMin, sharpenLabel(sharpenAll));

%% ------------------------------------------------------------------------
%% 3. CHAIN + SHARPEN + GATE
%% ------------------------------------------------------------------------
E = load(fullfile(cfgNom.dirs.energy, cfgNom.fname('energy', factor)));
Xk = E.X;  Uk = E.U;  sg = E.sigma;  tauf0 = E.tauf0;
rows = {};

for kr = 2:numel(rungs)                       % rung 1 IS the backbone
    T   = rungs(kr);
    cfg = minfuel_config(struct('thrustN', T));
    p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
    tf  = factor * tfmin(T);
    tag = sprintf('tulip_T%s_f%04d_eps%s', thrust_tag(T), round(1000*factor), ...
                  strrep(sprintf('%g', epsMin), '.', 'p'));
    matF = fullfile(outDir, [tag '.mat']);

    if isfile(matF) && ~rerun
        S = load(matF);  Xk = S.out.X(1:8,:);  Uk = S.out.U;
        fprintf('[%5.1f mN] cached -- skipping (rerun=true to redo)\n', T*1e3);
        rows{end+1} = S.row;
        continue
    end

    % --- energy continuation step (eps = 1) --------------------------------
    Xk(8,:) = Xk(8,:) * (tf / Xk(8,end));       % hold the factor: t_f,min ~ 1/T
    isLast  = (kr == numel(rungs));
    sched   = 1;
    if sharpenAll || isLast, sched = [1, cfg.schedSharpen(cfg.schedSharpen > epsMin), epsMin]; end

    okAll = true;  out = [];
    for e = sched
        o = struct('epsilon', e, 'tfTarget', tf, 'moonZone', -1, ...
                   'pSund', cfg.pSund, 'maxIter', 1500, 'warmTight', e < 1);
        out = casadi_energy_freetf(sg, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                                   Xk, Uk, tauf0, o);
        if ~(out.success && out.maxDefect < 1e-6)
            fprintf('[%5.1f mN] FAILED at eps=%.4g (defect %.2e) -- ladder stops\n', ...
                    T*1e3, e, out.maxDefect);
            okAll = false;  break
        end
        Xk = out.X(1:8,:);  Uk = out.U;
    end
    if ~okAll, break; end

    % --- gate --------------------------------------------------------------
    fp   = cr3bp_fingerprint(p, struct('tf', tf, 'factor', factor, 'epsMin', epsMin));
    gate = local_gate(out, doGate);
    dV   = p.c * log(1/out.mf) * p.lStar / p.tStar;
    row  = struct('T_mN', T*1e3, 'tf', tf, 'mf', out.mf, 'dV', dV, ...
                  'prop_kg', cfg.m0kg*(1-out.mf), 'switches', out.switches, ...
                  'defect', out.maxDefect, 'epsReached', sched(end), ...
                  'gatePass', gate.pass);
    rows{end+1} = row;

    sigma = sg;  %#ok<NASGU>
    save(matF, 'out', 'sigma', 'tauf0', 'rv0', 'rvf', 'T', 'fp', 'gate', 'row');
    fprintf('[%5.1f mN] eps=%g ok  defect=%.2e  m_f=%.6f  dV=%.4f  sw=%d  gate=%d\n', ...
            T*1e3, sched(end), out.maxDefect, out.mf, dV, out.switches, gate.pass);
end

%% ------------------------------------------------------------------------
%% 4. REPORT
%% ------------------------------------------------------------------------
fprintf('\n  T (mN)   t_f (ND)     m_f       dV (km/s)  prop (kg)  sw   defect    gate\n');
fprintf('  ---------------------------------------------------------------------------\n');
lines = {};
for k = 1:numel(rows)
    r = rows{k};
    s = sprintf('  %6.2f  %8.4f  %.6f   %7.4f    %7.4f  %3d  %.1e   %d', ...
                r.T_mN, r.tf, r.mf, r.dV, r.prop_kg, r.switches, r.defect, r.gatePass);
    fprintf('%s\n', s);  lines{end+1} = s;
end
sumF = fullfile(outDir, 'ladder_summary.txt');
fid = fopen(sumF, 'w');
fprintf(fid, 'tulip thrust ladder  %s\nfactor %.3f  epsMin %g  stepFrac %.3f\n\n', ...
        char(datetime('now','Format','yyyy-MM-dd HH:mm')), factor, epsMin, stepFrac);
fprintf(fid, '%s\n', lines{:});  fclose(fid);
fprintf('\n  summary -> %s\n', sumF);

% ---------------------------------------------------------------------------
function s = sharpenLabel(sharpenAll)
% SHARPENLABEL  Text for the header banner.
% INPUTS:  sharpenAll [logical]   OUTPUTS: s [char]
if sharpenAll, s = 'every rung'; else, s = 'final rung only'; end
end

% ---------------------------------------------------------------------------
function g = local_gate(out, doGate)
% LOCAL_GATE  The gate battery that is VALID for a free-time solution.
%
% Applied:
%   defect      - collocation residual, machine-tight required
%   boundSat    - no state/control box is saturated (a saturated box means the
%                 solution is shaped by an artificial bound, not the physics)
%   primer      - mean primer/thrust-direction misalignment, PMP first-order
%   lamMassEnd  - free-final-mass transversality residual
%
% NOT applied, and the reason matters:
%   certify_minfuel_pmp -- its arc-by-arc state/costate propagation assumes
%     dt/dtau = kappa. These solutions carry dt/dtau = c*kappa (the cScale slack
%     state), so the propagation would be wrong by the factor c and would return
%     a confidently misleading verdict. Needs a c-aware variant first.
%   run_foc_tulip -- it re-solves through casadi_minfuel_sundman with the
%     returnModel hook, at the campaign's NOMINAL thrust. It cannot gate an
%     off-nominal free-time rung.
%
% So a gatePass here is weaker than the nominal campaign's. It is stated rather
% than implied.
%
% INPUTS:  out - solver output struct; doGate - run the battery [logical]
% OUTPUTS: g   - struct with .pass and the individual measurements
g = struct('pass', false, 'defect', out.maxDefect, 'boundSatWorst', 'n/a', ...
           'primerDeg', NaN, 'lamMassEnd', NaN, 'notApplied', ...
           {{'certify_minfuel_pmp (assumes dt/dtau=kappa)', ...
             'run_foc_tulip (re-solves at nominal thrust)'}});
if ~doGate, return; end
% EXPECTED saturations, excluded from the verdict. `massHi` is the mass upper
% bound, and m(0) is pinned at the normalized maximum by the initial condition,
% so it is saturated at node 1 BY CONSTRUCTION -- every physically correct
% solution hits it. Failing on it made the gate report 0 on every healthy rung
% of the first ladder run (measured minSlack ~3e-12 on all five). A bound that
% every valid solution touches carries no information; the ones worth failing on
% are the artificial boxes (vBox, rBox, throttle) that would mean the answer is
% shaped by a limit rather than by the physics.
EXPECTED_SAT = {'massHi'};
if isfield(out, 'boundSat') && isstruct(out.boundSat)
    g.boundSatWorst = out.boundSat.worst;
    g.boundSatHit   = out.boundSat.hit && ~ismember(out.boundSat.worst, EXPECTED_SAT);
    g.boundSatExpected = ismember(out.boundSat.worst, EXPECTED_SAT);
else
    g.boundSatHit = false;  g.boundSatExpected = false;
end
if isfield(out, 'primerAlignDeg'), g.primerDeg  = out.primerAlignDeg; end
if isfield(out, 'lamMassEnd'),     g.lamMassEnd = out.lamMassEnd;     end
g.pass = out.maxDefect < 1e-6 && ~g.boundSatHit;
end
