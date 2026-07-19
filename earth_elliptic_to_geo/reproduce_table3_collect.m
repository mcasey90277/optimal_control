function tbl = reproduce_table3_collect(thrustList)
% REPRODUCE_TABLE3_COLLECT  Assemble + print the reproducer engine's
% "updated (best-found) Table 3" from already-computed REPRO_row_T*.mat
% files (reproduce_row.m's own output cache), one file per rung.
%
% Pure reader: never solves anything, never calls reproduce_row.m, no
% side effects beyond printing to stdout. For each requested thrust level,
% loads results/repro/REPRO_row_T<round(10*T)>.mat (reproduce_row.m's own
% cache filename convention -- see its `save(...)` line) and prints:
%   1. the fixed-width Table-3 block for that row (gergaud_row_str.m,
%      reused verbatim, not reimplemented);
%   2. a one-line comparison against the campaign FLOOR
%      (table3_certified.m): the delta in final mass, and a verdict of
%      "BEAT" (reproduced mass clearly above the floor), "matched" (within
%      a small numerical band of the floor), or -- defensively -- a
%      WARNING if the cached row is somehow BELOW the floor. This should
%      never happen for a row reproduce_row.m itself wrote (its own
%      verify_row.m call already throws on that breach before saving), but
%      this collector re-checks independently rather than trusting that
%      invariant blindly (e.g. a hand-edited or stale cache file).
% A rung whose REPRO_row_T*.mat file does not exist yet prints a MISSING
% line (NOT an error) and is OMITTED from the returned struct array -- this
% function only returns rows it actually found; it never fabricates a
% placeholder row for a missing rung.
%
% A requested thrust level with no registered campaign floor (e.g. the
% seeded-but-not-yet-attained 0.2/0.1 N deep rungs in table3_recipes.m) is
% still loaded and printed (row block only), with a NOTE line in place of
% the floor comparison, rather than erroring the whole collection.
%
% Also prints the R0-law spread: for every LOADED rung whose anchor was a
% REAL, independently-solved min-time anchor (anchor.anchorSource ~=
% 'R0law' -- i.e. 'solved' [coldB/chain strategies] or 'smallN_first'), the
% product T*tfmin (Gergaud's near-constant R0-law quantity; table3_recipes.m's
% own R0const=223.14) is printed per rung, plus its spread (min/max/range)
% across those rungs -- a sanity diagnostic that the independently
% re-solved anchors are mutually consistent with the R0 scaling law. Purely
% informational: not asserted/gated.
%
% INPUTS:
%   thrustList - vector of thrust levels [N] to collect, in the order to
%                print/return (e.g. [10 5 2.5 1 0.5])                    [1xK]
%
% OUTPUTS:
%   tbl - struct array of gergaud_row.m row structs, one entry per rung
%         that WAS found under results/repro/, in thrustList order, with
%         MISSING rungs omitted (numel(tbl) <= numel(thrustList); empty
%         struct array if none were found)                          [1xM struct]
%
% REFERENCES:
%   [1] .superpowers/sdd/task-4-brief.md (this function's spec).
%   [2] earth_elliptic_to_geo/reproduce_row.m (writes the REPRO_row_T*.mat
%       files this function reads: the results/repro/ path, the
%       REPRO_row_T<round(10*T)>.mat filename convention, and the
%       row/anchor/sol/rep saved-variable set).
%   [3] earth_elliptic_to_geo/gergaud_row_str.m / table3_certified.m /
%       table3_recipes.m (row formatting, campaign floor, R0const source).

if nargin < 1 || isempty(thrustList)
    error('reproduce_table3_collect:badInput', 'thrustList is required');
end

here     = fileparts(mfilename('fullpath'));
reproDir = fullfile(here, 'results', 'repro');

matchTolKg = 0.5;   % mirrors reproduce_row.m's own defaultTol(T).m_f_kg policy
                    % (0.5 kg numerical-noise slack below the campaign floor)

fprintf('\n=== Table 3 (updated, best-found) reproduction ===\n\n');

tbl      = struct([]);   % returned rows, found rungs only, in request order
anchorT  = [];            % thrust [N] of rungs with a REAL (non-R0law) anchor
anchorTf = [];            % matching tfmin [ND], same order as anchorT

for T = thrustList
    fname = fullfile(reproDir, sprintf('REPRO_row_T%d.mat', round(10*T)));
    if ~isfile(fname)
        fprintf('T=%-6g N: MISSING (%s not found -- run reproduce_row(%g) first)\n\n', ...
            T, fname, T);
        continue;
    end

    d    = load(fname, 'row', 'anchor');
    row  = d.row;
    anch = d.anchor;

    fprintf('%s', gergaud_row_str(row));

    try
        cert  = table3_certified(T);
        delta = row.m_f_kg - cert.m_f_kg;
        if delta > 1e-2
            fprintf(['  vs campaign floor: m_f=%.2f kg (floor %.2f kg) --> ' ...
                     '+%.2f kg BEAT | sw %d vs %d | revs %.3f vs %.3f\n'], ...
                row.m_f_kg, cert.m_f_kg, delta, row.switches, cert.switches, ...
                row.revs, cert.revs);
        elseif delta >= -matchTolKg
            fprintf(['  vs campaign floor: m_f=%.2f kg (floor %.2f kg) --> ' ...
                     'matched | sw %d vs %d | revs %.3f vs %.3f\n'], ...
                row.m_f_kg, cert.m_f_kg, row.switches, cert.switches, ...
                row.revs, cert.revs);
        else
            fprintf(['  *** WARNING *** vs campaign floor: m_f=%.2f kg is %.2f kg ' ...
                     'BELOW the certified floor %.2f kg (sw %d vs %d | revs %.3f vs ' ...
                     '%.3f) -- the reproducer should never regress below the floor; ' ...
                     'investigate this cached row\n'], ...
                row.m_f_kg, -delta, cert.m_f_kg, row.switches, cert.switches, ...
                row.revs, cert.revs);
        end
    catch
        fprintf(['  NOTE: no registered campaign floor for T=%g N (table3_certified ' ...
                 'has no entry -- likely a seeded/deep rung); printing row only, no ' ...
                 'floor comparison\n'], T);
    end
    fprintf('\n');

    if isempty(tbl)
        tbl = row;
    else
        tbl(end+1) = row; %#ok<AGROW>
    end

    if ~strcmp(anch.anchorSource, 'R0law')
        anchorT(end+1)  = T;          %#ok<AGROW>
        anchorTf(end+1) = anch.tfmin; %#ok<AGROW>
    end
end

fprintf('=== R0-law spread (T*tfmin, R0const=223.14) across real-anchor rungs ===\n');
if isempty(anchorT)
    fprintf('  (no loaded rungs had a real, independently-solved anchor)\n');
else
    prodTtf = anchorT .* anchorTf;
    for k = 1:numel(anchorT)
        fprintf('  T=%-6g N: tfmin=%.4f ND -> T*tfmin=%.3f\n', ...
            anchorT(k), anchorTf(k), prodTtf(k));
    end
    fprintf('  spread: min=%.3f max=%.3f range=%.3f (n=%d rungs)\n', ...
        min(prodTtf), max(prodTtf), max(prodTtf) - min(prodTtf), numel(prodTtf));
end
fprintf('\n');

end
