function S = sheet_from_arcs(arcs, opts)
%% Purpose:
%
%   Assemble the certified ARRIVAL-PHASE sheet (one departure phase, one
%   engine) from the grid crossings of any number of arclength_ms arcs.
%   Each crossing is a candidate root at an unwrapped level; the level is
%   folded back onto the grid (mod 1), duplicates reached by different arcs
%   are merged (same t_f within tolDup), EVERY candidate is certified with
%   certify_crossing and KEPT with its verdict, and S.TF(j) is the minimum
%   t_f over the certified candidates at grid point j -- NaN when none.
%   Nothing is discarded: a grid point can hold several local minima and
%   the failure reasons are part of the record.
%
%% Inputs:
%
%  arcs                     cell of struct          arclength_ms outputs (or
%                                                   anything with .crossings)
%  opts                     struct (optional)
%   .seeds struct array of already-certified candidates (certify_root
%   outputs with .sA) entered before the crossings -- the arcs' own start
%   points, which are not crossings of their own arc;
%   .sA [] an EXPLICIT list of arrival phases in [0,1), any spacing (the
%   sheet's columns, in this order); else the lattice .sA0 [0.0754] + (0:nA-1)/nA,
%   .nA [12];  .tolDup [1e-6] (days) and .tolZ [1e-6]
%   (relative, on z8) -- a duplicate must match in BOTH;
%   .certFn [@(p,sA) certify_crossing(p,sA,B,anc,copts)] -- default needs
%   .B, .anc (from arclength_arrival setup) and optional .copts;
%   .logFile ''
%
%% Outputs:
%
%  S                        struct                  .sA [1 x nA] .TF [1 x nA]
%                                                   days (NaN = none
%                                                   certified) .Z8 [8 x nA]
%                                                   .cand {1 x nA} struct
%                                                   arrays (certify_crossing
%                                                   outputs + .level .arc)
%                                                   .nCand .nCert .sA0 .nA
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
sA0 = d('sA0', 0.0754);  nA = d('nA', 12);  tolDup = d('tolDup', 1e-6);
% THE GRID IS A LIST. A lattice is the common case, but the columns may be
% any phases the caller names (run_phase_torus); every lookup below goes
% through the list, so nothing assumes a spacing.
if isfield(opts, 'sA') && ~isempty(opts.sA)
    sAlist = mod(opts.sA(:).', 1);  nA = numel(sAlist);  sA0 = sAlist(1);
else
    sAlist = mod(sA0 + (0:nA-1)/nA, 1);
end
assert(numel(unique(round(sAlist*1e9))) == nA, 'sheet_from_arcs: repeated arrival phases in the list');
tolZ = d('tolZ', 1e-6);
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
if isfield(opts, 'certFn') && ~isempty(opts.certFn)
    certFn = opts.certFn;
else
    assert(isfield(opts, 'B') && isfield(opts, 'anc'), 'sheet_from_arcs: pass .B and .anc or a .certFn');
    copts = d('copts', struct());
    certFn = @(p, sA) certify_crossing(p, sA, opts.B, opts.anc, copts);
end

S = struct('sA', sAlist, 'TF', nan(1, nA), 'Z8', nan(8, nA), ...
           'cand', {cell(1, nA)}, 'nCand', 0, 'nCert', 0, 'sA0', sA0, 'nA', nA);

% ---- pre-certified SEED points -----------------------------------------
% The arcs' own starting solutions are never crossings of their own arc, so
% a sheet built from crossings alone reports NaN at the very phases whose
% solutions launched it (measured: grid point 1 read "no certified
% candidate" while the 17.798 d anchor sat certified on disk). Seeds enter
% first, through the same merge and the same minimum.
seeds = d('seeds', struct([]));
for k = 1:numel(seeds)
    C = seeds(k);
    j = gridIndex(C.sA, sAlist);
    C.level = C.sA;  C.arc = 0;
    S.nCand = S.nCand + 1;
    if isempty(S.cand{j}), S.cand{j} = C; else, S.cand{j} = mergeStruct(S.cand{j}, C); end
    lg('  seed at j = %2d (sA %.4f): t_f %.6f d  %s', j, S.sA(j), C.tfDays, C.reason);
    if C.ok
        S.nCert = S.nCert + 1;
        if isnan(S.TF(j)) || C.tfDays < S.TF(j), S.TF(j) = C.tfDays;  S.Z8(:, j) = C.z(:); end
    end
end

for ia = 1:numel(arcs)
    cr = arcs{ia}.crossings;
    for ic = 1:numel(cr)
        c = cr(ic);
        j = gridIndex(c.level, sAlist);         % refuses a level off the list
        sA = S.sA(j);
        if ~c.converged
            C = struct('ok', false, 'reason', sprintf('crossing not converged (|R| = %.1e)', c.normR), ...
                       'z', nan(8,1), 'tfDays', NaN, 'sA', sA);
        else
            C = certFn(c.p, sA);
        end
        C.level = c.level;  C.arc = ia;
        % Merge with an existing candidate at this grid point only if it is
        % the SAME ROOT. Keying on final time alone would silently discard
        % a genuinely distinct extremal that happens to share a t_f -- the
        % grid point would lose a candidate and nothing would say so. The
        % costates settle it: same t_f AND same z8.
        dup = false;
        for k = 1:numel(S.cand{j})
            e = S.cand{j}(k);
            % A REFUSED candidate must never suppress a later SUCCESSFUL one:
            % one that failed a late gate still carries a finite t_f, so an
            % ok-blind merge could drop the certificate that fills the cell.
            if ~isfinite(C.tfDays) || abs(e.tfDays - C.tfDays) >= tolDup, continue, end
            if ~isequal(logical(e.ok), logical(C.ok)), continue, end
            if numel(e.z) == numel(C.z) && norm(e.z(:) - C.z(:)) <= tolZ*max(norm(C.z(:)), 1)
                dup = true;  break
            end
        end
        if dup
            lg('  arc %d crossing %d at j = %2d (sA %.4f): duplicate of an existing candidate (t_f %.6f d, same z8)', ia, ic, j, sA, C.tfDays);
            continue
        end
        S.nCand = S.nCand + 1;
        if isempty(S.cand{j}), S.cand{j} = C; else, S.cand{j} = mergeStruct(S.cand{j}, C); end
        lg('  arc %d crossing %d at j = %2d (sA %.4f): t_f %.6f d  %s', ia, ic, j, sA, C.tfDays, C.reason);
        if C.ok
            S.nCert = S.nCert + 1;
            if isnan(S.TF(j)) || C.tfDays < S.TF(j), S.TF(j) = C.tfDays;  S.Z8(:, j) = C.z(:); end
        end
    end
end
lg('sheet_from_arcs: %d candidates, %d certified, %d/%d grid points filled', ...
   S.nCand, S.nCert, nnz(isfinite(S.TF)), nA);
end

function j = gridIndex(level, sAlist)
% GRIDINDEX  Column of a (possibly unwrapped) phase level in the list:
% the entry within 1e-6 of it mod 1, or an error.  INPUTS: level; sAlist.
% OUTPUTS: j.
[dm, j] = min(abs(mod(sAlist - level + 0.5, 1) - 0.5));
assert(dm < 1e-6, 'phase %.6f is not on the grid', level);
end

function A = mergeStruct(A, C)
% MERGESTRUCT  Append C to struct array A, aligning fields.  INPUTS: A; C.
% OUTPUTS: A.
fa = fieldnames(A);  fc = fieldnames(C);
for k = 1:numel(fc), if ~isfield(A, fc{k}), [A.(fc{k})] = deal([]); end, end
for k = 1:numel(fa), if ~isfield(C, fa{k}), C.(fa{k}) = []; end, end
C = orderfields(C, A);
A(end+1) = C;
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function logmsg(f, s)
% LOGMSG  Print, and append to a log file when one is named.  INPUTS: f; s.
fprintf('%s\n', s);
if ~isempty(f), fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
