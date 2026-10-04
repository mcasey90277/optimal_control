function items = recertify_candidates(harvestMat, chunk, nChunk, outMat, opts)
%% Purpose:
%
%   The RE-CERTIFICATION share of the status-layer backfill (spec 6, steps 1
%   and 3): run the current certifier (certify_root: full stack, gates and
%   H6 after a conjugate failure) on one chunk of the work a harvest
%   (backfill_status_layer 'harvest') listed, saving after every item.
%
%     kind 'cand' -- a legacy candidate whose status is ambiguous (verdict
%                    0), re-seeded from its OWN junctions at its own phases;
%     kind 'stop' -- a rib's stall point, seeded from the rib's last
%                    certified point and solved in one step at the phase the
%                    rib stalled stepping to;
%     kind 'prim' -- a primary no source gave junctions, re-polished from its
%                    own z8 (seed_from_z8, K = 24).
%
%  ASSUMPTIONS / NOTES:
%
% • Item order: every 'cand' (needRecert order), then every 'stop', then
%   every 'prim' (unmatched order); this chunk takes items chunk:nChunk:end.
% • RESUME: an item already saved in outMat without an error is skipped; an
%   item saved WITH an error (the certifier threw) is retried and replaced.
% • moved = the result is not same_root as the item's own root ('cand': the
%   candidate's z, 'prim': the entry's z8). A moved primary is recorded and
%   never written (assemble refuses to use it); 'stop' is a new solve and
%   never moved.
% • The harvest key (backfill_status_layer 'key') is checked against
%   harvest.mat, stamped into outMat, and a resume refuses an outMat made
%   against another harvest (recertify_candidates:harvestKey).
% • Saves atomically: <outMat>.part then movefile. One log line per item in
%   <outMat>.log.
% • The fence: certify_root gets opts.pool (default capped_pool(1)); with no
%   pool it refuses unless opts.allowUnfenced. The seams .certifier .B
%   .seedFromZ8 exist for the test (no solve).
%
%% Inputs:
%
%  harvestMat               char                    harvest.mat (.cands
%                                                   .needRecert .stops
%                                                   .unmatched .unmatchedInfo
%                                                   .recordMat)
%  chunk, nChunk            double                  this chunk, of how many
%  outMat                   char                    this chunk's recert_<k>.mat
%  opts                     struct (optional)       .pool [capped_pool(1)]
%                                                   .allowUnfenced [false]
%                                                   .certifier [@certify_root]
%                                                   .B [] setup (default from
%                                                   the record) .seedFromZ8 []
%                                                   @(z8, rv0, K) -> seed
%                                                   .K [24] .copts [struct()]
%                                                   further certify_root options
%
%% Outputs:
%
%  items                    struct array            .kind .index .C (certify_root
%                                                   result) .moved .err ('' or
%                                                   the throw's message)
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 5, opts = struct(); end
H = load(harvestMat);
key = backfill_status_layer('key', H);
assert(isfield(H, 'harvestKey') && strcmp(H.harvestKey, key), 'recertify_candidates:harvestKey', ...
       '%s does not match its own stamped key: rebuild it', harvestMat);
certifier = fieldd(opts, 'certifier', @certify_root);
K = fieldd(opts, 'K', 24);
allowUnfenced = isfield(opts, 'allowUnfenced') && isequal(opts.allowUnfenced, true);
if isfield(opts, 'pool'), pool = opts.pool; else, pool = capped_pool(1); end
B = fieldd(opts, 'B', []);
if isempty(B)
    L = load(H.recordMat);  c = L.(char(fieldnames(L)));
    [B, ~] = arclength_arrival('setup', catalog_setup_request(c, c.sheets(1).sD_frac(1)));
end
seeder = fieldd(opts, 'seedFromZ8', @(z8, rv0, K_) seed_from_z8(z8, rv0, K_, B.Tnd, B.cnd, B.mu));
copts = fieldd(opts, 'copts', struct());
copts.pool = pool;  copts.allowUnfenced = allowUnfenced;

% ---- the work list and this chunk's share -----------------------------------
stops = struct([]);  if isfield(H, 'stops'), stops = H.stops; end
info = struct([]);   if isfield(H, 'unmatchedInfo'), info = H.unmatchedInfo; end
nC = numel(H.needRecert);  nS = numel(stops);  nP = numel(info);
kinds = [repmat({'cand'}, 1, nC), repmat({'stop'}, 1, nS), repmat({'prim'}, 1, nP)];
idx = [H.needRecert(:)', 1:nS, [info.k]];
pos = [1:nC, 1:nS, 1:nP];                    % position within its kind's list
mine = chunk:nChunk:numel(kinds);

% ---- resume ----------------------------------------------------------------------
items = struct('kind', {}, 'index', {}, 'C', {}, 'moved', {}, 'err', {});
if isfile(outMat)
    R = load(outMat);
    if ~(isfield(R, 'harvestKey') && strcmp(R.harvestKey, key))
        error('recertify_candidates:harvestKey', ['%s was made against another harvest (key mismatch): its indices ' ...
              'do not address these candidates; move it away'], outMat);
    end
    items = R.items;
end
harvestKey = key;
logF = [outMat '.log'];

for q = mine
    kind = kinds{q};  index = idx(q);
    at = find(strcmp({items.kind}, kind) & [items.index] == index, 1);
    if ~isempty(at) && isempty(items(at).err), continue, end
    switch kind
        case 'cand'
            Cc = H.cands(index);  sD = Cc.sD;  sA = Cc.sA;  zRef = Cc.z(:);  Yref = Cc.Y;
        case 'stop'
            s = stops(pos(q));  sD = s.sDto;  sA = s.sA;  zRef = s.z(:);  Yref = s.Y;
        case 'prim'
            p = info(pos(q));  sD = p.sD;  sA = p.sA;  zRef = p.z8(:);  Yref = [];
    end
    t0 = tic;  err = '';  C = struct();  moved = false;
    try
        rv0 = B.stateD(sD);  rvf = B.stateA(sA);  rv0 = rv0(1:6);  rvf = rvf(1:6);
        if isempty(Yref)
            seed = seeder(zRef, rv0, K);
        else                                 % the junctions it carries, node 1 forced
            Kc = size(Yref, 2);
            seed = struct('tf', zRef(8), 'tGrid', linspace(0, zRef(8), Kc + 1), 'Y', [Yref, Yref(:, end)]);
            seed.Y(1:7, 1) = [rv0(:); 1];  seed.Y(8:14, 1) = zRef(1:7);
        end
        co = copts;  co.sD = mod(sD, 1);  co.sA = mod(sA, 1);
        C = certifier(seed, rv0, rvf, B, co);
        moved = ~strcmp(kind, 'stop') && ~same_root(C.z, zRef);
    catch ME
        err = ME.message;
    end
    it = struct('kind', kind, 'index', index, 'C', C, 'moved', moved, 'err', err);
    if isempty(at), items(end+1) = it; else, items(at) = it; end
    part = [outMat '.part'];
    save(part, 'items', 'harvestKey', '-mat');  movefile(part, outMat, 'f');
    fid = fopen(logF, 'a');
    fprintf(fid, '%s %s %d at (%.6f, %.6f): %s | %.0f s\n', char(datetime('now', 'Format', 'HH:mm:ss')), kind, index, ...
            sD, sA, verdict(C, moved, err), toc(t0));
    fclose(fid);
end
end

% ---------------------------------------------------------------------------
function s = verdict(C, moved, err)
% VERDICT  One-line outcome.  INPUTS: C; moved; err.  OUTPUTS: s.
if ~isempty(err), s = ['THREW: ' err];  return, end
st = NaN;  if isfield(C, 'status'), st = C.status; end
r = '';  if isfield(C, 'reason'), r = C.reason; end
s = sprintf('status %g%s -- %s', st, tern(moved, ' MOVED', ''), r);
end

function s = tern(c, a, b)
% TERN  Ternary.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
