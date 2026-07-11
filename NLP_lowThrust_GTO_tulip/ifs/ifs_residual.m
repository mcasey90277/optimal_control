function [R, J] = ifs_residual(Z, prob)
% IFS_RESIDUAL  Switch-structured multiple-shooting residual for IFS.
%
% Integrates each hard-throttle arc and assembles the square residual
%   R = [continuity_1..k (16 each); terminal (8); switch_1..k (1 each)],
% continuity_a = arcEndpoint_a - N_a, terminal = BC on the last arc's endpoint
% (rendezvous or fixedState), switch_i = S(N_i). Jacobian (Task 3) is a per-arc
% complex-step sparse block matrix.
%
% INPUTS:
%   Z    - unknown vector [(8+17k)x1] (see IFS_PACK)
%   prob - problem struct (see plan Shared data layout)
% OUTPUTS:
%   R - residual [(8+17k)x1]
%   J - Jacobian dR/dZ [(8+17k)x(8+17k)] sparse (Task 3; [] here)
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
k = prob.k;
[lam0, N, tau] = ifs_unpack(Z, k);

% arc start states and tau-spans
[starts, spans] = ifs_arcs(lam0, N, tau, prob);   % starts{a}[16], spans{a}[1x2]
e = cell(1, k+1);
for a = 1:k+1
    e{a} = ifs_int_arc(starts{a}, spans{a}, prob.uArc(a), prob);
end

Rcont = zeros(16, k);
for a = 1:k
    Rcont(:, a) = e{a} - N(:, a);
end
Rterm = ifs_termres(e{k+1}, prob);
Rsw = zeros(k, 1);
for ii = 1:k
    Rsw(ii) = ifs_S(N(:, ii), prob.c);
end
R = [Rcont(:); Rterm; Rsw];

if nargout > 1
    J = [];   % filled in Task 3
end
end

% ---- helpers shared with the Jacobian (Task 3) --------------------------
function [starts, spans] = ifs_arcs(lam0, N, tau, prob)
% Build per-arc start state and tau-span from the unknowns.
k = prob.k;  starts = cell(1,k+1);  spans = cell(1,k+1);
starts{1} = [prob.rv0(:); prob.m0; prob.t0; lam0];
spans{1}  = [prob.tau0, tau(1)];            % arc-1 tau-span starts at prob.tau0
for a = 2:k
    starts{a} = N(:, a-1);
    spans{a}  = [tau(a-1), tau(a)];
end
starts{k+1} = N(:, k);
spans{k+1}  = [tau(k), prob.tauf];
end

function eE = ifs_int_arc(startY, span, uArc, prob)
% Integrate one arc; returns the endpoint [16x1] (may be complex under CS).
[~, Y] = ode113(@(s,y) ifs_eom(s, y, prob.Tmax, prob.c, prob.muStar, ...
                              prob.pSund, uArc), span, startY, prob.odeOpts);
eE = Y(end, :).';
end

function S = ifs_S(Y, c)
% Switching function at a node state.
S = 1 - sqrt(sum(Y(12:14).^2))*c/Y(7) - Y(15);
end

function Rt = ifs_termres(eE, prob)
% Terminal residual [8x1] on the last arc's endpoint.
if strcmp(prob.termMode, 'rendezvous')
    Rt = [eE(1:6) - prob.rvf(:); eE(15); eE(8) - prob.tf];   % rv, lamM=0, t=tf
else   % fixedState
    Rt = eE(1:8) - prob.termTarget(:);
end
end
