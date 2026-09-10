function ok = test_certify_caps()
%% Purpose:
%
%   Tests the HARD WALL-CLOCK FENCE on certify_root's external calls.
%
%   Why this test exists: on 2026-09-09 the sheet assembly and a departure
%   rib both hung -- 12 and 16 hours at 100% CPU, no progress, inside ONE
%   function evaluation. certify_root passed `wallSec` to ms_tfmin and I
%   took that for protection. It is not: an in-process wall check fires
%   only BETWEEN solver iterations, and a CR3BP segment propagation whose
%   iterate parks a junction near a primary can grind for hours inside a
%   single evaluation. That is documented in run_capped's own header, added
%   after the same thing happened twice on the 0.5 N campaign. The only
%   real fence is a parfeval cap that cancels the worker.
%
%   A capped call must FAIL LOUDLY with a named reason. Silently accepting
%   a candidate whose witness or gates could not be run would put an
%   unverified entry in the catalog, which is the one thing the gate stack
%   exists to prevent.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

pool = capped_pool();
ok = chk(ok, ~isempty(pool) && isvalid(pool), sprintf('pool available (%d workers)', pool.NumWorkers));

[B, anc] = arclength_arrival('setup');

% (1) generous caps: no regression, the anchor still certifies
C = certify_crossing(anc.p, anc.sA, B, anc, struct('pool', pool));
ok = chk(ok, C.ok && abs(C.tfDays - 17.7976) < 2e-3, ...
         sprintf('fenced but generous: anchor certifies at %.4f d (%s)', C.tfDays, C.reason));

% (2) each external call, capped to nothing, must REFUSE and name itself
caps = {'capPolishSec', 'polish'; 'capFlySec', 'flight'; ...
        'capWitnessSec', 'witness'; 'capGatesSec', 'gates'};
for k = 1:size(caps, 1)
    o = struct('pool', pool);  o.(caps{k,1}) = 1e-4;
    Ck = certify_crossing(anc.p, anc.sA, B, anc, o);
    hit = ~Ck.ok && contains(lower(Ck.reason), 'cap');
    ok = chk(ok, hit, sprintf('%s cap refuses with a named reason: %s', caps{k,2}, Ck.reason));
end

if ok, fprintf('TEST_CERTIFY_CAPS: ALL PASS\n'); else, fprintf('TEST_CERTIFY_CAPS: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
