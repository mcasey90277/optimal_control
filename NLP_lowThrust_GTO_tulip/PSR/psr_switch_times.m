function [sigSw, arcU, sigGrid] = psr_switch_times(solFile)
% PSR_SWITCH_TIMES  Bang-bang switch locations (in the Sundman mesh variable).
%
% Returns the throttle switch points of a bang-bang solution as sub-node
% crossings of s = 0.5 on the solution's own sigma mesh, plus the per-arc
% throttle pattern. Shared by the movie title and psr_switch_hessian.
%
% INPUTS:
%   solFile - solution .mat (or struct) in seed layout: out.U [4xnN], sigma
% OUTPUTS:
%   sigSw   - [1xk] switch locations in the sigma variable (sub-node, linear)
%   arcU    - [1x(k+1)] throttle in {0,1} on each arc (between consecutive
%             switches; arc 1 = [sigma(1), sigSw(1)], etc.)
%   sigGrid - [1xnN] the sigma mesh (row)
%
% REFERENCES: PSR/psr_second_order.m, PSR/run_psr.m stage 6.

if ischar(solFile) || isstring(solFile), S = load(solFile); else, S = solFile; end
s   = S.out.U(4, :);
sig = S.sigma(:).';
cr  = find(diff(sign(s - 0.5)) ~= 0);        % node index before each crossing
sigSw = zeros(1, numel(cr));
for q = 1:numel(cr)
    k = cr(q);
    sigSw(q) = sig(k) + (0.5 - s(k))*(sig(k+1) - sig(k))/(s(k+1) - s(k));
end
edges = [sig(1), sigSw, sig(end)];
arcU  = zeros(1, numel(sigSw) + 1);
for a = 1:numel(arcU)
    mid = 0.5*(edges(a) + edges(a+1));
    arcU(a) = double(interp1(sig, s, mid, 'linear') > 0.5);
end
sigGrid = sig;
end
