function ok = test_dro_tulip_seed()
%% Purpose:
%
%   Tests dro_tulip_seed -- "which certified solution seeds this phase
%   pair, and is it the right engine?", the lookup transfer_study section 4
%   used to spell out inline.
%
%   Checks:
%     1. a file-backed pair and a catalog-only pair both return the seed
%        seed_from_entry builds for that entry, BITWISE;
%     2. info names the entry that supplied it;
%     3. an off-grid phase pair is refused as :noSeed;
%     4. a DIFFERENT operating point is refused as :operatingPoint --
%        separately, so the two failures are distinguishable (the inline
%        version had one assert for both and could not say which);
%     5. the operating point is checked in ALL six fields, not just the two
%        phases: each one perturbed alone must refuse.
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
if ~isfile(fullfile(ind, 'results', 'costate_catalog_dro_tulip_70mN.mat'))
    fprintf('  SKIP  no 70 mN catalog on disk\n');  return
end

ndp = nd_propulsion(0.070, 900, 150);   mu = 0.012150585609624;
B = arclength_arrival('setup');
op = struct('tau', 1.0, 'Np', 7, 'pm', -1, 'thrustN', 0.070, 'ispS', 900, 'm0kg', 150, ...
            'Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', mu);
lib = dro_tulip_library(ind, struct('includeCatalog', true));

pairs = {0, 0.0754, 'file-backed anchor'; 1/12, 0.0754 + 2/12, 'catalog-only'};
for k = 1:2
    sD = pairs{k,1};  sA = pairs{k,2};
    o = op;  o.rv0 = B.stateD(sD);
    [seed, info] = dro_tulip_seed(sD, sA, o);
    j = find(abs([lib.sD] - mod(sD,1)) < 1e-6 & abs([lib.sA] - mod(sA,1)) < 1e-6, 1);
    ref = seed_from_entry(lib(j), o.rv0, o);
    ok = chk(ok, isequal(seed.Y, ref.Y) && isequal(seed.tGrid, ref.tGrid) && seed.tf == ref.tf, ...
             sprintf('%s pair: the seed is seed_from_entry''s, bitwise', pairs{k,3}));
    ok = chk(ok, strcmp(info.src, lib(j).src) && abs(info.sA - lib(j).sA) < 1e-12, ...
             sprintf('   and info names the entry (%s)', info.src));
end

o = op;  o.rv0 = B.stateD(0);
ok = chk(ok, refuses(@() dro_tulip_seed(0, 0.1, o), 'dro_tulip_seed:noSeed'), ...
         'an off-grid phase pair is refused as :noSeed');

perturb = {'tau', 2.0; 'Np', 5; 'pm', 1; 'thrustN', 0.015; 'ispS', 1710; 'm0kg', 15};
allSix = true;
for k = 1:size(perturb, 1)
    bad = o;  bad.(perturb{k,1}) = perturb{k,2};
    allSix = allSix && refuses(@() dro_tulip_seed(0, 0.0754, bad), 'dro_tulip_seed:operatingPoint');
end
ok = chk(ok, allSix, 'each of the six operating-point fields refuses on its own');

if ok, fprintf('TEST_DRO_TULIP_SEED: ALL PASS\n'); else, fprintf('TEST_DRO_TULIP_SEED: FAIL\n'); end
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
