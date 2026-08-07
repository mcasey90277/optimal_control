function val = vizOption(opts,name,val)
%% Purpose:
%
%  Return one option's value when the caller supplied it and the default
%  otherwise. Every function in coorbital.viz resolves its options through
%  this, so the option block of each reads as one unbranched line per option,
%  the same shape the entry scripts' override blocks have.
%
%  An ABSENT field and a field set to [] are treated alike, both meaning "use
%  the default". That equivalence is what lets a caller write
%  struct('Parent',[]) to mean the same thing as omitting Parent, which in
%  turn is what makes a programmatically assembled option struct easy to build
%  -- a caller can set every field unconditionally and leave the ones it does
%  not care about empty.
%
%% Inputs:
%
%  opts             Struct                      Caller options
%
%  name             Char [1 x n]                Option name
%
%  val              Any                         Default value
%
%% Outputs:
%
%  val              Any                         opts.(name) when that field is
%                                               present and non-empty,
%                                               otherwise the default
%
%% Notes:
%
%  A private function, so it carries no self-demo: nothing outside
%  +coorbital/+viz can call it, and inside the package it is reached by plain
%  name rather than by namespace.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    if isfield(opts,name) && ~isempty(opts.(name))
               val = opts.(name);
    end
end
