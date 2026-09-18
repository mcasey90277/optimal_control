function ok = test_certify_schema()
% TEST_CERTIFY_SCHEMA  Every certificate has ONE field set, whichever gate it
% stopped at. Producers append certificates into struct arrays (the sheet's
% seeds, a rib's points), and MATLAB refuses to assign a struct into an array
% with different fields -- so a field added only on the success path aborts
% the sheet at its first mixed outcome (found in review 2026-09-18, after the
% H2/H3 bounds were added that way).
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
pool = capped_pool(2);
lib = dro_tulip_library();  E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup', struct('physicsOnly', true));
rv0 = B.stateD(E.sD);  rvf = B.stateA(E.sA);
seed = seed_from_z8(E.z(:), rv0(1:6), 24, B.Tnd, B.cnd, B.mu);
base = struct('sA', E.sA, 'sD', E.sD, 'pool', pool);

Cgood = certify_root(seed, rv0(1:6), rvf(1:6), B, base);
Cearly = certify_root(seed, rv0(1:6), rvf(1:6), B, setfield(base, 'gateKm', 1e-9));      %#ok<SFLD> stops at the flight gate
Chyp = certify_root(seed, rv0(1:6), rvf(1:6), B, setfield(base, 'hypFloor', 1e9));       %#ok<SFLD> stops at H2
ok = chk(ok, Cgood.ok && ~Cearly.ok && ~Chyp.ok, sprintf('one pass, one early refusal (%s), one H2 refusal', Cearly.reason));
ok = chk(ok, isequal(sort(fieldnames(Cgood)), sort(fieldnames(Cearly))) && isequal(sort(fieldnames(Cgood)), sort(fieldnames(Chyp))), ...
         'all three certificates carry the same fields');
try
    arr = Cgood;  arr(end+1) = Cearly;  arr(end+1) = Chyp; %#ok<NASGU>
    ok = chk(ok, true, 'they append into one struct array, as the sheet and the ribs do');
catch ME
    ok = chk(ok, false, ['they do NOT append into one struct array: ' ME.message]);
end
ok = chk(ok, isfinite(Cgood.minLamVBound) && isfinite(Cgood.minQmtBound) && isnan(Cearly.minLamVBound), ...
         'the H2/H3 bounds are numbers on a pass and NaN where the gate was never reached');
if ok, fprintf('test_certify_schema: ALL PASS\n'); else, fprintf('test_certify_schema: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
