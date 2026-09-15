function [ok, msg, info] = rib_validate(f, spec)
%% Purpose:
%
%   Decide whether a rib file is a PUBLISHABLE result for a SPECIFIC unit.
%   The first version checked that a variable named `problem` existed and
%   counted whatever sat in .pts; it would have accepted another column's
%   rib, an uncertified point, or an empty identity (Astra pass 3, 3.4).
%   This one checks structure, identity, coordinates and certification
%   flags against the spec the campaign hands it, and reports coverage.
%
%   Coverage is reported, not required: a walker that stalls before the
%   last grid point returns a real, terminal, SHORTER column, and the
%   library ships what is certified -- but the caller must see that.
%
%% Inputs:
%
%  f                        char                    rib .mat path
%  spec                     struct (optional)       .nPts points requested
%                                                   per rib; .col expected
%                                                   R.j; .sA expected R.sA;
%                                                   .nD departure lattice;
%                                                   .problem identity struct
%                                                   to match (thrustN ispS
%                                                   m0kg tauDRO NpTulip
%                                                   pmTulip sD). A numeric
%                                                   spec is taken as .nPts.
%
%% Outputs:
%
%  ok                       logical                 publishable
%  msg                      char                    why not, or ''
%  info                     struct                  .nPts .stop .complete
%                                                   .sD (certified phases)
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, spec = struct(); end
if isnumeric(spec), spec = struct('nPts', spec); end
d = @(fld, v) fieldd(spec, fld, v);
nPts = d('nPts', NaN);
ok = false;  msg = '';  info = struct('nPts', 0, 'stop', '', 'complete', false, 'sD', []);
try
    if ~isfile(f), msg = sprintf('%s does not exist', f);  return, end
    w = whos('-file', f);
    if ~any(strcmp({w.name}, 'R')), msg = sprintf('%s holds no rib struct R', f);  return, end
    L = load(f, '-mat');                 % a worker's temporary is named .part: without -mat, load reads it as ASCII
    R = L.R;
    if ~isstruct(R) || isempty(R), msg = 'rib struct R is empty';  return, end
    if numel(R) ~= 1, msg = sprintf('rib file holds %d ribs; a unit is exactly one', numel(R));  return, end
    for req = {'j', 'sA', 'pts', 'stop'}
        if ~isfield(R, req{1}), msg = sprintf('rib struct has no .%s', req{1});  return, end
    end
    if isfinite(d('col', NaN)) && R.j ~= spec.col
        msg = sprintf('rib is column %d, this unit is column %d', R.j, spec.col);  return
    end
    if isfinite(d('sA', NaN)) && abs(mod(R.sA - spec.sA + 0.5, 1) - 0.5) > 1e-9
        msg = sprintf('rib arrival phase %.6f is not the unit''s %.6f', R.sA, spec.sA);  return
    end
    info.stop = char(string(R.stop));
    P = R.pts;
    if isempty(P), msg = sprintf('rib has no certified points (stop: %s)', info.stop);  return, end
    if ~isstruct(P) || ~all(isfield(P, {'ok', 'sD', 'z'}))
        msg = 'rib points lack .ok/.sD/.z';  return
    end
    if ~all([P.ok]), msg = sprintf('%d of %d rib points are not certified (.ok false)', nnz(~[P.ok]), numel(P));  return, end
    sD = [P.sD];
    if any(~isfinite(sD)), msg = 'a rib point has a non-finite departure phase';  return, end
    sDlist = d('sD', []);
    if ~isempty(sDlist)                      % an explicit departure list
        sDlist = mod(sDlist(:).', 1);
        off = arrayfun(@(v) min(abs(mod(sDlist - v + 0.5, 1) - 0.5)), sD);
        if any(off > 1e-6), msg = sprintf('%d rib point(s) are off the %d-phase departure list', nnz(off > 1e-6), numel(sDlist));  return, end
        if numel(unique(round(mod(sD, 1)*1e9))) ~= numel(sD), msg = 'rib has duplicate departure phases';  return, end
    end
    nD = d('nD', NaN);
    if isfinite(nD) && isempty(sDlist)
        off = abs(mod(sD*nD + 0.5, 1) - 0.5);
        if any(off > 1e-6*nD), msg = sprintf('%d rib point(s) are off the %d-point departure lattice', nnz(off > 1e-6*nD), nD);  return, end
        if numel(unique(round(mod(sD, 1)*nD))) ~= numel(sD), msg = 'rib has duplicate departure phases';  return, end
    end
    if ~any(strcmp({w.name}, 'problem')), msg = 'rib carries no problem identity';  return, end
    Pr = L.problem;
    if ~isstruct(Pr) || isempty(fieldnames(Pr)), msg = 'rib problem identity is empty';  return, end
    want = d('problem', []);
    if isstruct(want)
        for fn = {'thrustN', 'ispS', 'm0kg', 'tauDRO', 'NpTulip', 'pmTulip', 'sD'}
            if ~isfield(want, fn{1}), continue, end
            if ~isfield(Pr, fn{1}), msg = sprintf('rib identity has no %s', fn{1});  return, end
            if abs(Pr.(fn{1}) - want.(fn{1})) > 1e-9*max(1, abs(want.(fn{1})))
                msg = sprintf('rib was certified at %s = %g, this campaign is %g', fn{1}, Pr.(fn{1}), want.(fn{1}));  return
            end
        end
    end
    info.nPts = numel(P);  info.sD = sD;
    info.complete = strcmp(info.stop, 'complete') && (isnan(nPts) || info.nPts == nPts);
    ok = true;
catch ME
    ok = false;  msg = sprintf('rib file could not be validated: %s', ME.message);
end
end

function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present and non-empty, else v0.  INPUTS: s; f; v0.
% OUTPUTS: v.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
