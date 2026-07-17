function r = mod(a, b)
% MOD  Local override of MATLAB's built-in mod(), extended to accept CasADi
% MX/SX symbolic arguments.
%
% Needed because lt_mee_rhs.m (Task 1, verified, not modified here) calls
% mod(L, 2*pi), and Task 3's collocation makes L = pi + sigma*DeltaL a
% symbolic MX expression once DeltaL is a decision variable. CasADi's MX
% class does not overload MATLAB's built-in mod() (only its own .fmod()/
% .remainder(), which use C-fmod truncated/sign-of-dividend semantics, not
% MATLAB's sign-of-divisor mod() semantics) -- so mod(MX, double) errors with
% "Invalid data type" as a plain builtin call. This shim leaves ALL numeric
% behavior byte-identical (delegates straight to builtin('mod',...)) and only
% takes the symbolic branch for casadi.MX/SX inputs, via the textbook
% identity that IS the built-in's own documented definition:
% mod(a,m) = a - m.*floor(a./m). floor() is supported by CasADi MX/SX.
% Placed in this folder (not edited into lt_mee_rhs.m, which stays the
% verified Task 1 source of truth) -- MATLAB's name resolution shadows the
% built-in for every caller while this folder is on the path, including
% lt_mee_rhs.m, so no upstream file changes.
%
% INPUTS:
%   a - dividend [any size; double, or casadi.MX/SX]
%   b - divisor  [any size; double, or casadi.MX/SX]
%
% OUTPUTS:
%   r - a mod b, MATLAB mod() semantics [size of a/b per implicit expansion]
%
% REFERENCES:
%   [1] MATLAB documentation for mod(): "mod(a,m) = a - m.*floor(a./m))"
if isa(a, 'casadi.MX') || isa(a, 'casadi.SX') || isa(b, 'casadi.MX') || isa(b, 'casadi.SX')
    r = a - b.*floor(a./b);
else
    r = builtin('mod', a, b);
end
end
