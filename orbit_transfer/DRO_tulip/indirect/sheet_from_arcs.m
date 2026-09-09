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
%   .sA0 [0.0754] .nA [12] grid;  .tolDup [1e-6] (days) duplicate merge;
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
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
if isfield(opts, 'certFn') && ~isempty(opts.certFn)
    certFn = opts.certFn;
else
    assert(isfield(opts, 'B') && isfield(opts, 'anc'), 'sheet_from_arcs: pass .B and .anc or a .certFn');
    copts = d('copts', struct());
    certFn = @(p, sA) certify_crossing(p, sA, opts.B, opts.anc, copts);
end

S = struct('sA', mod(sA0 + (0:nA-1)/nA, 1), 'TF', nan(1, nA), 'Z8', nan(8, nA), ...
           'cand', {cell(1, nA)}, 'nCand', 0, 'nCert', 0, 'sA0', sA0, 'nA', nA);

for ia = 1:numel(arcs)
    cr = arcs{ia}.crossings;
    for ic = 1:numel(cr)
        c = cr(ic);
        j = mod(round((c.level - sA0)*nA), nA) + 1;         % grid index
        assert(abs(mod(c.level - sA0, 1)*nA - round(mod(c.level - sA0, 1)*nA)) < 1e-6 || ...
               abs(mod(c.level - sA0, 1)*nA - nA) < 1e-6, 'crossing level %.6f is not on the grid', c.level);
        sA = S.sA(j);
        if ~c.converged
            C = struct('ok', false, 'reason', sprintf('crossing not converged (|R| = %.1e)', c.normR), ...
                       'z', nan(8,1), 'tfDays', NaN, 'sA', sA);
        else
            C = certFn(c.p, sA);
        end
        C.level = c.level;  C.arc = ia;
        % merge with an existing candidate at this grid point (same root)
        dup = false;
        for k = 1:numel(S.cand{j})
            if isfinite(C.tfDays) && abs(S.cand{j}(k).tfDays - C.tfDays) < tolDup, dup = true; break, end
        end
        if dup
            lg('  arc %d crossing %d at j = %2d (sA %.4f): duplicate of an existing candidate (t_f %.6f d)', ia, ic, j, sA, C.tfDays);
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
