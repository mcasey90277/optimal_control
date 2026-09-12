function ok = test_seed_from_entry()
%% Purpose:
%
%   Tests seed_from_entry -- one library entry becomes one multiple-shooting
%   seed. The same construction existed THREE times: inline in
%   transfer_study section 4, as run_dro_tulip's private seed_of, and inside
%   build_arrival_sheet's seeding loop.
%
%   Two routes, because a library entry comes in two shapes:
%     FILE-BACKED  it carries junction states already -- reuse them, append
%                  the final column the solver wants, and pin column 1 to
%                  the ACTUAL departure state and the entry's costates;
%     CATALOG      it carries z8 only -- fly it once and cut (seed_from_z8).
%
%   Checks, all against the inline constructions they replace:
%     1. the file-backed route is BITWISE what the three sites built;
%     2. column 1 is the actual departure state and the entry's costates,
%        not whatever the stored trajectory began with;
%     3. the catalog route is BITWISE seed_from_z8;
%     4. the seed contract holds: tGrid spans [0, tf], K+1 columns;
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

lib = dro_tulip_library(ind, struct('includeCatalog', true));
ndp = nd_propulsion(0.070, 900, 150);   mu = 0.012150585609624;
B = arclength_arrival('setup');
kFile = find(arrayfun(@(e) ~isempty(e.Y), lib), 1);
kCat  = find(arrayfun(@(e)  isempty(e.Y), lib), 1);
P = lib(kFile);   Q = lib(kCat);
rv0P = B.stateD(P.sD);   rv0Q = B.stateD(Q.sD);

% 1 + 2. the file-backed route, against the inline construction
K = size(P.Y, 2);
old = struct('tf', P.z(8), 'tGrid', linspace(0, P.z(8), K+1), 'Y', [P.Y, P.Y(:,end)]);
old.Y(1:7,1) = [rv0P(1:6); 1];   old.Y(8:14,1) = P.z(1:7);
new = seed_from_entry(P, rv0P, struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', mu));
ok = chk(ok, isequal(new.Y, old.Y) && isequal(new.tGrid, old.tGrid) && new.tf == old.tf, ...
         sprintf('file-backed entry [%s] reproduces the inline seed BITWISE', P.src));
ok = chk(ok, isequal(new.Y(1:7,1), [rv0P(1:6); 1]) && isequal(new.Y(8:14,1), P.z(1:7)), ...
         'column 1 is the actual departure state and the entry''s costates');

% 3. the catalog route, against seed_from_z8
oldQ = seed_from_z8(Q.z, rv0Q(1:6), 24, ndp.Tnd, ndp.cnd, mu);
newQ = seed_from_entry(Q, rv0Q, struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', mu, 'K', 24));
ok = chk(ok, isequal(newQ.Y, oldQ.Y) && isequal(newQ.tGrid, oldQ.tGrid), ...
         sprintf('catalog entry (%.4f, %.4f) reproduces seed_from_z8 BITWISE', Q.sD, Q.sA));

% 4. the seed contract
ok = chk(ok, size(new.Y, 2) == numel(new.tGrid) && abs(new.tGrid(end) - new.tf) < 1e-15 && new.tGrid(1) == 0, ...
         'the seed contract holds: K+1 columns, grid spanning [0, tf]');

% 5. refusals
bad = Q;  bad.z = Q.z(1:7);
ok = chk(ok, refuses(@() seed_from_entry(bad, rv0Q, struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', mu)), ...
         'seed_from_entry:z8'), 'an entry whose z is not 8 long is refused');
ok = chk(ok, refuses(@() seed_from_entry(Q, rv0Q, struct()), 'seed_from_entry:physics'), ...
         'a catalog entry with no thrust/exhaust/mu to fly it is refused');

if ok, fprintf('TEST_SEED_FROM_ENTRY: ALL PASS\n'); else, fprintf('TEST_SEED_FROM_ENTRY: FAIL\n'); end
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
