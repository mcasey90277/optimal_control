function best = pilot_rung_20mN()
% PILOT_RUNG_20MN  Ladder-prep validation: one warm-chained 20 mN fuel rung
% on the ELFO free-t_f (freetf) engine.
%
% Chains the f1.200 GTO->ELFO min-energy seed to the new 20 mN thrust via
% chain_rung_seed_elfo (trivial pass-through + fingerprint; cScale decouples
% the clock so no resample is needed), GATES an explicit energy re-clean at
% the new thrust (casadi_energy_freetf, eps=1, t_f pinned to the seed's own
% t_f, cBox rung-scaled per spec sec 4), then sharpens eps 1->0 through
% gen_elfo_minfuel's hardened homotopy. PASS = the sharpen returns (it
% internally asserts full certification before saving) with a clean
% boundSat and fp recorded. Artifacts under _T20mN tags; certified caches
% (energy_elfo_f1200.mat) untouched. Both the re-clean gate and the sharpen
% call are wrapped so an honest failure prints a BLOCKED line and returns
% instead of crashing (house honesty rule; mirrors sundman_homotopy's soft
% .certified flag on the tulip side). (2026-07-21 ladder-prep T6; spec sec 6.)
%
% OUTPUTS: best - struct(.certified, .stage, .ipoptStatus, .maxDefect,
%          .switches, .boundSat, .fp) describing the pilot's final state
% REFERENCES: [1] spec 2026-07-21-ladder-prep-design.md sec 6.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();

cfgNom = minfuel_config();                       % nominal (25 mN) reference
cfg    = minfuel_config(struct('thrustN', 0.020));
p20    = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
tag    = thrust_tag(cfg.thrustN);                % '_T20mN'

seedFile = fullfile(here, 'results', 'energy_elfo_f1200.mat');
S = load(seedFile);
[S, fpChain] = chain_rung_seed_elfo(S, p20, struct('pilot','20mN'));

% --- gated energy re-clean at the new thrust (eps=1, t_f pinned) -----------
tfTarget = S.X(8,end);
Tfac  = cfgNom.thrustN / cfg.thrustN;             % 0.025/0.020 = 1.25
cBoxD = [0.15*min(1,Tfac), 6*max(1,Tfac)];        % [0.15 7.5]
oClean = struct('moonZone',S.moonZone, 'muGain',1, 'tfTarget',tfTarget, ...
    'epsilon',1, 'pSund',S.pSund, 'qSund',S.qSund, 'tfCapMult',6, ...
    'cBox',cBoxD, 'maxIter',3000, 'warmTight',true);
rClean = casadi_energy_freetf(S.sigma, S.rv0, S.rvf, p20.Tmax, p20.c, ...
                               p20.muStar, S.X, S.U, S.tauf0, oClean);
cleanOk = strcmp(rClean.ipoptStatus,'Solve_Succeeded') && rClean.maxDefect < 1e-6;
if ~cleanOk
    fprintf('\nPILOT 20mN ELFO: BLOCKED at energy re-clean (status=%s defect=%.2g)\n', ...
            rClean.ipoptStatus, rClean.maxDefect);
    best = struct('certified',false, 'stage','energy re-clean', ...
                   'ipoptStatus',rClean.ipoptStatus, 'maxDefect',rClean.maxDefect, ...
                   'switches',rClean.switches, 'boundSat',rClean.boundSat, 'fp',[]);
    return
end

% --- save the re-cleaned 20 mN seed WITH fp (_T20mN-tagged) ----------------
fp = cr3bp_fingerprint(p20, struct('tf', tfTarget, 'chainedFrom', fpChain.chainedFrom));
[~, baseName] = fileparts(seedFile);
seedT20File = fullfile(here, 'results', sprintf('%s%s.mat', baseName, tag));
X = rClean.X; U = rClean.U; sigma = S.sigma; rv0 = S.rv0; rvf = S.rvf;
tauf0 = S.tauf0;  moonZone = S.moonZone;  pSund = S.pSund;  qSund = S.qSund;
tf = X(8,end);
save(seedT20File, 'X','U','sigma','rv0','rvf','tauf0','moonZone','pSund','qSund','tf','fp');

% --- eps 1->0 sharpen on the re-cleaned seed (gen_elfo_minfuel's own gates) --
outFile = fullfile(here, 'results', sprintf('pilot_minfuel%s.mat', tag));
try
    gen_elfo_minfuel(struct('seedFile',seedT20File, 'target','ELFO', ...
        'outFile',outFile, 'maxIter',3000, 'thrustN',cfg.thrustN));
catch sharpErr
    fprintf('\nPILOT 20mN ELFO: BLOCKED at sharpen (%s: %s)\n', ...
            sharpErr.identifier, sharpErr.message);
    best = struct('certified',false, 'stage','sharpen', 'ipoptStatus',sharpErr.message, ...
                   'maxDefect',NaN, 'switches',NaN, 'boundSat',[], 'fp',fp);
    return
end

% gen_elfo_minfuel asserts full certification internally before saving, so a
% normal return here means the saved artifact IS certified.
R = load(outFile);
sat = 'n/a'; if isfield(R.out,'boundSat'), sat = R.out.boundSat.worst; end
nSw = sum(abs(diff(R.U(4,:) > 0.5)));
best = struct('certified',true, 'stage','sharpen', 'ipoptStatus',R.out.ipoptStatus, ...
              'maxDefect',R.out.maxDefect, 'switches',nSw, 'boundSat',R.out.boundSat, ...
              'fp',R.fp);
fprintf('\nPILOT 20mN ELFO: certified=%d status=%s defect=%.2g sw=%d boundSatWorst=%s\n', ...
        best.certified, best.ipoptStatus, best.maxDefect, best.switches, sat);
end
