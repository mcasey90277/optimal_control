function [ok, msg, info] = rib_validate(f, nPts)
%% Purpose:
%
%   Decide whether a rib file is a PUBLISHABLE unit result: it exists,
%   loads, holds a rib struct with certified points, and carries a problem
%   identity. Coverage is reported, not required: a walker that stalls
%   before the last grid point returns a real, terminal, SHORTER column
%   (rib_from_crossing's .stop says where and why), and the library ships
%   what is certified -- but the caller must be able to see that the column
%   is short rather than count its file as full coverage (Astra pass 2, F).
%
%% Inputs:
%
%  f                        char                    rib .mat path
%  nPts                     double (optional)       points requested per rib
%
%% Outputs:
%
%  ok                       logical                 publishable
%  msg                      char                    why not, or ''
%  info                     struct                  .nPts .stop .complete
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, nPts = NaN; end
ok = false;  msg = '';  info = struct('nPts', 0, 'stop', '', 'complete', false);
if ~isfile(f), msg = sprintf('%s does not exist', f);  return, end
try
    w = whos('-file', f);
catch ME
    msg = sprintf('%s does not load: %s', f, ME.message);  return
end
if ~any(strcmp({w.name}, 'R')), msg = sprintf('%s holds no rib struct R', f);  return, end
L = load(f, 'R');
if isempty(L.R) || ~isfield(L.R, 'pts'), msg = 'rib struct R is empty or has no .pts';  return, end
info.nPts = sum(arrayfun(@(r) numel(r.pts), L.R));
info.stop = strjoin(arrayfun(@(r) char(string(r.stop)), L.R, 'UniformOutput', false), '; ');
info.complete = all(arrayfun(@(r) strcmp(r.stop, 'complete'), L.R)) && (isnan(nPts) || info.nPts == nPts*numel(L.R));
if info.nPts == 0, msg = sprintf('rib has no certified points (stop: %s)', info.stop);  return, end
if ~any(strcmp({w.name}, 'problem')), msg = 'rib carries no problem identity';  return, end
ok = true;
end
