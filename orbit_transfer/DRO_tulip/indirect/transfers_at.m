function T = transfers_at(c, sD, sA, opts)
%% Purpose:
%
%   Every transfer the catalog records at one (departure phase, arrival
%   phase) pair, for a mission designer: the cell's primary plus all
%   alternatives within a circular phase tolerance, ranked by flight time,
%   each with its optimality status, reason, t_f [days], Delta-V,
%   propellant, z8, junctions and source.
%
%  ASSUMPTIONS / NOTES:
%
% • Single-sheet catalogs only (a named error otherwise).
% • A catalog WITHOUT the status layer returns its primary with status NaN
%   and statusName 'not annotated', and no alternatives.
% • Phases are matched circularly (s = 1 - 1e-9 finds s = 0).
% • Delta-V and propellant follow the catalog's own derive formulas:
%   Tmax = (rungs_N(1)/m0)*tStar^2/(lStar*1000), mf = 1 - Tmax*tf/c_nd,
%   dV = c_nd*log(1/mf)*lStar/tStar.
% • Nothing recorded -> an EMPTY struct array with the same fields.
%
%% Inputs:
%
%  c                        struct                  single-sheet costate catalog
%  sD, sA                   scalar                  departure / arrival phase
%                                                   fractions [0,1)
%  opts.tol                 scalar                  circular phase tolerance
%                                                   [1e-6]
%
%% Outputs:
%
%  T                        struct array            sorted by tfDays; fields
%                                                   kind sD sA status
%                                                   statusName reason tfDays
%                                                   dvKms propellantKg z8
%                                                   junctions source
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'tol'), opts.tol = 1e-6; end
if ~isscalar(c.sheets)
    error('transfers_at:multiSheet', 'transfers_at supports single-sheet catalogs only (got %d sheets).', numel(c.sheets));
end
K = status_key();
near = @(a, b) abs(mod(a - b + 0.5, 1) - 0.5) <= opts.tol;
sh = c.sheets(1);
rows = {};

% Primary: the sheet cell at (sD, sA)
iD = find(near(sh.sD_frac(:), sD), 1);
iA = find(near(sh.sA_frac(:), sA), 1);
if ~isempty(iD) && ~isempty(iA) && sh.has_solution(iD, iA)
    k = sh.entry_index(iD, iA);
    if isfield(sh, 'status') && ~isempty(sh.status)
        st = double(sh.status(iD, iA));  rs = sh.status_reason{k};  jn = sh.junctions{k};
    else
        st = NaN;  rs = '';  jn = [];
    end
    rows{end+1} = mkrow(c, K, 'primary', sh.sD_frac(iD), sh.sA_frac(iA), sh.z8(:, k), st, rs, jn, 'sheet');
end

% Alternatives
if isfield(c, 'alternatives') && ~isempty(c.alternatives)
    for ka = 1:numel(c.alternatives)
        a = c.alternatives(ka);
        if near(a.sD, sD) && near(a.sA, sA)
            rows{end+1} = mkrow(c, K, 'alternative', a.sD, a.sA, a.z8, double(a.status), a.status_reason, a.junctions, a.source);
        end
    end
end

if isempty(rows)
    T = struct('kind', {}, 'sD', {}, 'sA', {}, 'status', {}, 'statusName', {}, 'reason', {}, ...
               'tfDays', {}, 'dvKms', {}, 'propellantKg', {}, 'z8', {}, 'junctions', {}, 'source', {});
    return
end
T = [rows{:}];
[~, ord] = sort([T.tfDays]);
T = T(ord);
end

function r = mkrow(c, K, kind, sD, sA, z8, st, reason, junctions, source)
% MKROW  Build one output row.  INPUTS: catalog c, status key K, kind, phases,
% z8 [8x1], status, reason, junctions, source.  OUTPUTS: r [struct].
tf   = z8(8);
Tmax = (c.rungs_N(1)/c.thruster.m0_kg)*c.constants.tStar_s^2/(c.constants.lStar_km*1000);
mf   = 1 - Tmax*tf/c.thruster.c_nd;
r.kind = kind;  r.sD = sD;  r.sA = sA;  r.status = st;
if isnan(st)
    r.statusName = 'not annotated';
else
    r.statusName = K.names{K.codes == st};
end
r.reason = reason;
r.tfDays = tf*c.constants.tStar_s/86400;
r.dvKms = c.thruster.c_nd*log(1/mf)*c.constants.lStar_km/c.constants.tStar_s;
r.propellantKg = (1 - mf)*c.thruster.m0_kg;
r.z8 = z8;  r.junctions = junctions;  r.source = source;
end
