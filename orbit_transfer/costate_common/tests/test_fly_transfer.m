function ok = test_fly_transfer()
%% Purpose:
%
%   Tests fly_transfer -- fly converged costates ONCE, validate that flight,
%   and attach what every consumer then recomputes: final mass fraction,
%   propellant, and Delta-V from the rocket equation. transfer_study section
%   5, certify_root, audit_phase_catalog and the witness in
%   verify_with_pumpkyn each spelled this out.
%
%   Checks:
%     1. the flight is BITWISE what the inline block built (same propagator
%        call, same struct fields);
%     2. mass, propellant and Delta-V match the inline arithmetic exactly;
%     3. Delta-V agrees with catalog_schema's named `deltav_from_mf`
%        derivation on the real catalog -- two homes, one rocket equation;
%     4. an INADMISSIBLE flight comes back reported, not thrown: the caller
%        decides what a bad flight means (the script asserts, the certifier
%        returns a reason);
%     5. malformed input is refused by name.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
ind = fullfile(fileparts(here), 'DRO_tulip', 'indirect');
addpath(here, ind);
if ~isfile(fullfile(ind, 'results', 'mintime_70mN_anchor.mat'))
    fprintf('  SKIP  no 70 mN anchor on disk\n');  return
end

L = load(fullfile(ind, 'results', 'mintime_70mN_anchor.mat'));  z8 = L.z(:);
ndp = nd_propulsion(0.070, 900, 150);
lStar = 389703.264829278;  tStar = 382981.289129055;  mu = 0.012150585609624;
B = arclength_arrival('setup');
rv0 = B.stateD(0);  rvf = B.stateA(0.0754);
phys = struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', mu, ...
              'lStar', lStar, 'tStar', tStar, 'm0kg', 150);

% the inline block transfer_study section 5 used to carry
[tu, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6); 1; z8(1:7)], ndp.Tnd, ndp.cnd, mu);
mfRef = Y(end,7);
dvRef = ndp.cnd*log(1/mfRef)*lStar/tStar;
VFref = validate_flight(tu, Y, z8(8), ndp.Tnd, ndp.cnd, mu, lStar);

F = fly_transfer(z8, rv0(1:6), rvf(1:6), phys);
ok = chk(ok, isequal(F.t, tu) && isequal(F.Y, Y), 'the flight is bitwise the inline propagation');
ok = chk(ok, isequal(F.mf, mfRef) && isequal(F.dvKms, dvRef) && isequal(F.propellantKg, (1-mfRef)*150), ...
         sprintf('mass %.6f, dV %.4f km/s and propellant %.2f kg match the inline arithmetic', ...
                 F.mf, F.dvKms, F.propellantKg));
ok = chk(ok, F.admissibility.ok == VFref.ok && strcmp(F.admissibility.reason, VFref.reason), ...
         sprintf('admissibility carried through: %s', F.admissibility.reason));
ok = chk(ok, all(isfield(F, {'t','Y','z8','rv0','rvf','Tnd','cnd','muStar','lStar','tStar','m0kg','nSamples'})), ...
         'the flight carries the physics and endpoints downstream needs');

% 3. one rocket equation, two homes
catMat = fullfile(ind, 'results', 'costate_catalog_dro_tulip_70mN.mat');
if isfile(catMat)
    Lc = load(catMat);  fn = fieldnames(Lc);  cat_ = Lc.(fn{1});
    dvSchema = catalog_schema('derive', cat_, 'deltav_from_mf', struct('mf_frac', F.mf));
    ok = chk(ok, abs(dvSchema - F.dvKms) <= 4*eps(F.dvKms), ...
             sprintf('dV agrees with catalog_schema''s derivation (%.10f vs %.10f km/s)', dvSchema, F.dvKms));
end

% 4. a bad flight is REPORTED, not thrown
% scaled costates leave the trajectory admissible (the min-time direction is
% invariant under a positive scaling), so that would be a check that cannot
% fail. Fly far past mass depletion instead: the all-burn law m = 1 - T t/c
% goes negative, which is exactly what the validator exists to catch.
bad = z8;  bad(8) = z8(8)*5;
Fb = fly_transfer(bad, rv0(1:6), rvf(1:6), phys);
ok = chk(ok, isstruct(Fb) && isfield(Fb, 'admissibility') && ~Fb.admissibility.ok, ...
         sprintf('an inadmissible flight is REPORTED, not thrown: %s', Fb.admissibility.reason));

% 5. refusals
ok = chk(ok, refuses(@() fly_transfer(z8(1:7), rv0(1:6), rvf(1:6), phys), 'fly_transfer:z8'), ...
         'a z8 that is not 8 long is refused');
ok = chk(ok, refuses(@() fly_transfer(z8, rv0(1:6), rvf(1:6), rmfield(phys, 'cnd')), 'fly_transfer:physics'), ...
         'missing physics is refused');

if ok, fprintf('TEST_FLY_TRANSFER: ALL PASS\n'); else, fprintf('TEST_FLY_TRANSFER: FAIL\n'); end
end

function r = refuses(fh, id)
% REFUSES  True when fh throws the named identifier.  INPUTS: fh; id.
% OUTPUTS: r [logical].
try
    fh();  r = false;
catch ME
    r = strcmp(ME.identifier, id);
    if ~r, fprintf('        (threw %s: %s)\n', ME.identifier, ME.message); end
end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
