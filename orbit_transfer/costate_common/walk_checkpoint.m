function varargout = walk_checkpoint(action, f, varargin)
%% Purpose:
%
%   A resumable CHECKPOINT for a sequential walk (a rib column): the state
%   after the last accepted point, saved atomically, so that an attempt
%   killed nine hours into a column is continued by the next attempt from
%   the last accepted point instead of from zero (Astra pass 2/3: "a
%   reclaimed column restarts from zero" -- the declared loss budget was a
%   whole column per kill).
%
%   THE CHECKPOINT BELONGS TO A UNIT, NOT AN ATTEMPT: it lives at
%   <output>.ckpt, is written only by the unit's owner (the worker holding
%   its lock), and carries the unit's identity (column, arrival phase,
%   lattice, direction, targets, problem). A resumer VALIDATES it against
%   what it is about to walk and refuses a mismatch by name. The worker
%   removes it after the unit is published.
%
%% Inputs:
%
%  action                   char                    'save' | 'load' | 'clear'
%  f                        char                    checkpoint path
%  varargin                 save:  S (the state struct; .identity required)
%                           load:  expect (identity struct to match)
%
%% Outputs:
%
%  save  -> none
%  load  -> [S, why]: S = [] and why = reason when absent or mismatched
%  clear -> none
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

switch lower(action)
    case 'save'
        S = varargin{1};
        assert(isfield(S, 'identity'), 'walk_checkpoint:identity', 'a checkpoint must carry .identity');
        S.version = 1;  S.saved = char(datetime('now'));
        tmp = sprintf('%s.%s.part', f, char(java.util.UUID.randomUUID()));
        save(tmp, '-struct', 'S');
        publish_atomic(tmp, f);
        varargout = {};

    case 'load'
        expect = varargin{1};
        if ~isfile(f), varargout = {[], 'no checkpoint'};  return, end
        try
            C = load(f, '-mat');            % the extension is .ckpt: without -mat, load reads it as ASCII
        catch ME
            varargout = {[], sprintf('checkpoint unreadable: %s', ME.message)};  return
        end
        if ~isfield(C, 'version') || C.version ~= 1 || ~isfield(C, 'identity')
            varargout = {[], 'checkpoint has no version/identity'};  return
        end
        [same, what] = sameIdentity(C.identity, expect);
        if ~same, varargout = {[], sprintf('checkpoint is for another walk (%s)', what)};  return, end
        varargout = {C, ''};

    case 'clear'
        if isfile(f), delete(f); end
        varargout = {};

    otherwise
        error('walk_checkpoint:action', 'unknown action "%s"', action);
end
end

% ------------------------------------------------------------------------
function [same, what] = sameIdentity(a, b)
% SAMEIDENTITY  Field-by-field equality of two identity structs, numeric
% fields to 1e-9 relative, arrays by value.  INPUTS: a; b.  OUTPUTS: same;
% what (the first differing field).
same = false;  what = '';
fa = fieldnames(a);  fb = fieldnames(b);
if ~isempty(setxor(fa, fb)), what = sprintf('fields differ: %s', strjoin(setxor(fa, fb), ' '));  return, end
for k = 1:numel(fb)
    x = a.(fb{k});  y = b.(fb{k});
    if isstruct(x) && isstruct(y)
        [s2, w2] = sameIdentity(x, y);
        if ~s2, what = sprintf('%s.%s', fb{k}, w2);  return, end
    elseif isnumeric(x) && isnumeric(y)
        if ~isequal(size(x), size(y)) || any(abs(x(:) - y(:)) > 1e-9*max(1, abs(y(:)))), what = fb{k};  return, end
    elseif ~isequal(x, y)
        what = fb{k};  return
    end
end
same = true;
end
