function ver = verify_cr3bp_pmp(opts)
% VERIFY_CR3BP_PMP  Campaign driver: lunar-aware PMP verification for a
% certified CR3BP (elliptic->GEO WITH lunar gravity) front-door min-fuel
% artifact.
%
% Loads a run_cr3bp_geo.m front-door product (results/cr3bp_T<N>N_phi<phi>_
% fuel.mat: fields X [7x(N+1)], U [4x(N+1)], defectDuals [7xN], sigma
% [(N+1)x1], L, dL, fp incl. thrustN/m0kg/ispS/muM/DM/nM/phi0/gain), rebuilds
% the (par, par.pert, out, sigma) triple verify_pmp_mee.m needs, and runs
% the SAME 2-body verifier (earth_elliptic_to_geo/direct/verify/
% verify_pmp_mee.m) -- now lunar-aware per its 2026-07-23 amendment (task B):
% mee_primer_switch.m subtracts the zero-throttle ballistic/lunar bracket out
% of its B(X)/pel extraction and its S-formula G0 term before forming the
% primer vector and switching function, so the primer/sign gates are an
% HONEST check under lunar gravity, not the pre-amendment corrupted numbers
% solve_cr3bp_minfuel.m's CAVEAT recorded (certification there rested on the
% four NLP metrics alone, explicitly NOT on primer/switching agreement).
%
% This driver is a thin reproduction script, no new verification logic (house
% convention, mirrors earth_elliptic_to_geo/direct/verify/run_verify_pmp_mee.m):
% it only assembles inputs for, and prints the output of, verify_pmp_mee.m
% (+ switch_structure.m, a purely-primal dual-free cross-check, and
% hamiltonian_const_check.m for the time-costate-constancy diagnostic -- see
% the HAMILTONIAN CAVEAT below for why that one is reported, not gated).
%
% PASS/FAIL is the suite's OWN threshold (verify_pmp_mee.m's ver.pass:
% primerMedianDeg<1 deg AND overallSignPct>=99%) -- this driver does not
% invent or relax any gate.
%
% HAMILTONIAN CAVEAT: hamiltonian_const_check.m's constancy argument assumes
% the dynamics/cost are AUTONOMOUS in the time state t (Pontryagin's
% dlambda_t/dL = -dH/dt = 0). Under lunar gravity (par.pert active) the Moon's
% position is ang = pert.nM*t + pert.phi0 -- an EXPLICIT function of t -- so
% the CR3BP problem is genuinely non-autonomous and lambda_t is NOT expected
% to be a first integral. This driver still runs the check and prints its
% verdict for visibility, but a NOT-CONSTANT verdict here is the CORRECT
% physical answer, not a defect; it is reported, never gated.
%
% INPUTS:
%   opts - (optional) struct, all fields optional:
%     .thrustN - max thrust level [N], must match an existing front-door
%                artifact's filename tag (default 10)             [scalar]
%     .phi0    - lunar phase at t=0 [rad], filename tag (default 0) [scalar]
%     .matFile - explicit artifact path, overrides thrustN/phi0 filename
%                construction if given (default '' -> auto-built)   [char]
%
% OUTPUTS:
%   ver - verify_pmp_mee.m's output struct (primer/switching gate numbers +
%         per-node arrays), PLUS:
%     .matFile          - artifact path actually loaded                [char]
%     .fp               - the artifact's fingerprint struct           [struct]
%     .switchStructure  - switch_structure.m output (primal cross-check) [struct]
%     .hamiltonianCheck - hamiltonian_const_check.m output (see CAVEAT) [struct]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/direct/verify/verify_pmp_mee.m,
%       mee_primer_switch.m (the lunar-aware verifier this drives; see the
%       latter's header "LUNAR-AWARE AMENDMENT" for the full derivation).
%   [2] earth_elliptic_to_geo/direct/verify/run_verify_pmp_mee.m (2-body
%       reproduction-driver pattern this file mirrors).
%   [3] earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m (produces this
%       driver's input artifact; field names/units match verbatim).
%   [4] earth_elliptic_to_geo_CR3BP/direct/solve_cr3bp_minfuel.m (recorded
%       the CAVEAT this task closes -- see its header comment).
%   [5] earth_elliptic_to_geo_CR3BP/TODO.md ("CR3BP-aware primer + PSR" item;
%       this file closes the primer half).
if nargin < 1 || isempty(opts), opts = struct(); end
setup_paths();   % adds this folder + the 2-body campaign's core/lib/verify

d = @(f, v) optdef(opts, f, v);
thrustN = d('thrustN', 10);
phi0    = d('phi0', 0);
matFile = d('matFile', '');

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if isempty(matFile)
    matFile = fullfile(resDir, sprintf('cr3bp_T%sN_phi%s_fuel.mat', ...
        num_tag(thrustN), num_tag(phi0)));
end
assert(exist(matFile, 'file') == 2, 'verify_cr3bp_pmp:missingArtifact', ...
    'front-door artifact not found: %s -- run run_cr3bp_geo.m first', matFile);

S  = load(matFile);
fp = S.fp;
assert(isfield(fp, 'gain') && fp.gain > 0, 'verify_cr3bp_pmp:notLunar', ...
    ['artifact %s has gain=%.4g -- this driver targets LUNAR-AWARE ' ...
     '(gain>0) verification; a gain=0 (pure 2-body) artifact is already ' ...
     'covered by the earth_elliptic_to_geo campaign''s own verify suite'], ...
    matFile, fp.gain);

% --- rebuild par + par.pert from the artifact's own fingerprint (fp) --------
% (fp.muM/DM/nM/phi0/gain are the EXACT lunar_params.m outputs the solve used
% -- read directly, not recomputed via lunar_params, so this driver can never
% silently diverge from the solved physics via a units/rounding mismatch.)
par = kepler_lt_params(fp.thrustN, fp.m0kg, fp.ispS);
par.pert = struct('muM', fp.muM, 'DM', fp.DM, 'nM', fp.nM, ...
                   'phi0', fp.phi0, 'gain', fp.gain);

out   = struct('X', S.X, 'U', S.U, 'dL', S.dL, 'lamDef', S.defectDuals);
sigma = S.sigma;

fprintf('\n=== verify_cr3bp_pmp: %s ===\n', matFile);
fprintf('    T=%g N, m0=%g kg, Isp=%g s, gain=%.4f (lunar mass scale), phi0=%.4f rad\n', ...
    fp.thrustN, fp.m0kg, fp.ispS, fp.gain, fp.phi0);

% --- primer/switching-function verification (the pert-aware suite) ---------
ver = verify_pmp_mee(out, par, sigma, struct('eps', 0));

% --- primal (dual-free) switch-structure cross-check ------------------------
sw = switch_structure(S.X, S.U, S.dL, S.sigma);
fprintf(['[switch_structure] nSwitch=%d | revs=%.3f | duty=%.4f | ' ...
         'nodesPerRev=%.2f\n'], sw.nSwitch, sw.revs, sw.duty, sw.nodesPerRev);

% --- Hamiltonian (time-costate constancy) diagnostic; see HAMILTONIAN CAVEAT
resWrap = struct('fuel', struct('lamDef', S.defectDuals, 'dL', S.dL), 'sigma', S.sigma);
chk = hamiltonian_const_check(resWrap);

passStr = 'FAIL';  if ver.pass, passStr = 'PASS'; end
fprintf(['\nverify_cr3bp_pmp: PRIMER/SWITCH GATE = %s  ' ...
         '(primerMedianDeg=%.3f [<1 required], overallSignPct=%.2f%% [>=99 required])\n'], ...
    passStr, ver.primerMedianDeg, ver.overallSignPct);
fprintf(['verify_cr3bp_pmp: Hamiltonian verdict = %s (INFORMATIONAL under lunar ' ...
         'gravity -- the Moon''s motion makes the problem non-autonomous in t, ' ...
         'so lambda_t is NOT expected to be a first integral here; see HAMILTONIAN ' ...
         'CAVEAT in this file''s header)\n'], chk.verdict);

ver.matFile          = matFile;
ver.fp               = fp;
ver.switchStructure  = sw;
ver.hamiltonianCheck = chk;
end

% ---------------------------------------------------------------------------
function s = num_tag(v)
% NUM_TAG  Filename-safe numeric tag: integers -> plain digits ('10'),
% non-integers -> decimal point replaced by 'p' ('0.5' -> '0p5'), negative
% sign replaced by 'm'. Replicated locally from run_cr3bp_geo.m /
% solve_cr3bp_minfuel.m's helper of the same name (campaign convention: kept
% local so this file has no cross-file helper dependency).
if abs(v - round(v)) < 1e-9
    s = sprintf('%d', round(v));
else
    s = strrep(sprintf('%g', v), '.', 'p');
end
s = strrep(s, '-', 'm');
end
