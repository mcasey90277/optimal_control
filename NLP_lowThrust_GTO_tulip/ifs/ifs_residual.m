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
%   J - Jacobian dR/dZ [(8+17k)x(8+17k)] sparse, per-arc complex-step
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
k = prob.k;
[lam0, N, gblock] = ifs_unpack(Z, k);

% arc start states and tau-spans (gblock -> tau via ifs_taus inside ifs_arcs)
[starts, spans] = ifs_arcs(lam0, N, gblock, prob);   % starts{a}[16], spans{a}[1x2]
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
    nR = 8 + 17*k;  hCS = 1e-20;
    contRow = @(a) (a-1)*16 + (1:16);    % continuity block a
    termRow = 16*k + (1:8);
    swRow   = @(ii) 16*k + 8 + ii;
    nodeIdx = @(ii) 8 + (ii-1)*16 + (1:16);
    tauIdx  = @(ii) 8 + 16*k + ii;
    Jt = [];  Ji = [];  Jj = [];         % triplet accumulators (filled via pushblock)

    % --- per-arc complex-step blocks (integrated sensitivities) ----------
    % NB: ifs_eom is autonomous (tau unused), so a complex tau-span endpoint is
    % a valid analytic continuation -> CS through the integration limit is exact.
    for a = 1:k+1
        % which Z-unknowns arc a's integration depends on:
        if a == 1
            zdep = (1:8).';                       % lam0 (state fixed)
        else
            zdep = nodeIdx(a-1).';                % start node N_{a-1}
        end
        for jg = 1:min(a, k)                            % stick-breaking: arc a depends on g_1..g_min(a,k)
            zdep = [zdep; tauIdx(jg)];
        end
        % residual rows this arc's endpoint feeds:
        if a <= k
            rr = contRow(a);  applyG = @(de) de;                  % Rcont_a = e_a - N_a
        else
            rr = termRow;     applyG = @(de) ifs_dterm(de, prob); % Rterm = g(e_{k+1})
        end
        for c = zdep.'
            Zp = Z;  sc = max(1, abs(Zp(c)));  Zp(c) = Zp(c) + 1i*hCS*sc;
            [sP, spP] = ifs_arcs_one(Zp, a, prob);
            eP = ifs_int_arc(sP, spP, prob.uArc(a), prob);
            de = imag(eP)/(hCS*sc);
            [Ji, Jj, Jt] = pushblock(Ji, Jj, Jt, rr, c, applyG(de));
        end
    end
    % --- analytic blocks -------------------------------------------------
    for a = 1:k                          % Rcont_a = e_a - N_a  -> -I on N_a
        [Ji, Jj, Jt] = pushblock(Ji, Jj, Jt, contRow(a), nodeIdx(a), -eye(16));
    end
    for ii = 1:k                         % Rsw_i = S(N_i) -> dS/dN_i
        dS = ifs_dS(N(:, ii), prob.c);   % [1x16]
        [Ji, Jj, Jt] = pushblock(Ji, Jj, Jt, swRow(ii), nodeIdx(ii), dS);
    end
    J = sparse(Ji, Jj, Jt, nR, nR);
end
end

% ---- helpers shared with the Jacobian (Task 3) --------------------------
function [starts, spans] = ifs_arcs(lam0, N, gblock, prob)
% Build per-arc start state and tau-span from the unknowns. gblock is the
% unconstrained gap-param block, mapped to bounded switch times via IFS_TAUS.
k = prob.k;  starts = cell(1,k+1);  spans = cell(1,k+1);
tau = ifs_taus(gblock, prob.tau0, prob.tauf);
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

function [Ji, Jj, Jt] = pushblock(Ji, Jj, Jt, rrows, ccols, B)
% PUSHBLOCK  Accumulate a dense block B into triplet lists for every
% (row,col) pair, for later assembly via sparse(Ji,Jj,Jt,...).
%
% INPUTS:
%   Ji, Jj, Jt - running triplet lists (row idx, col idx, value) [vectors]
%   rrows      - residual row indices this block occupies [vector]
%   ccols      - unknown column indices this block occupies [vector]
%   B          - dense sensitivity block, size numel(rrows) x numel(ccols)
% OUTPUTS:
%   Ji, Jj, Jt - triplet lists with B's entries appended [vectors]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
[RR, CC] = ndgrid(rrows(:), ccols(:));
Ji = [Ji; RR(:)];  Jj = [Jj; CC(:)];  Jt = [Jt; B(:)];
end

function [startY, span] = ifs_arcs_one(Z, a, prob)
% IFS_ARCS_ONE  Rebuild just arc a's start state and tau-span from a
% (possibly complex, under complex-step) unknown vector Z.
%
% INPUTS:
%   Z    - unknown vector [(8+17k)x1], possibly complex (CS perturbation)
%   a    - arc index [scalar]
%   prob - problem struct (see plan Shared data layout)
% OUTPUTS:
%   startY - arc a's start state [16x1], possibly complex
%   span   - arc a's tau-span [1x2], possibly complex
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
k = prob.k;  [lam0, N, gblock] = ifs_unpack(Z, k);
[starts, spans] = ifs_arcs(lam0, N, gblock, prob);
startY = starts{a};  span = spans{a};
end

function dg = ifs_dterm(de, prob)
% IFS_DTERM  Map an endpoint sensitivity de[16] to terminal-residual
% sensitivity [8x1], per the terminal-mode selection matrix.
%
% INPUTS:
%   de   - endpoint state sensitivity [16x1]
%   prob - problem struct (see plan Shared data layout)
% OUTPUTS:
%   dg   - terminal-residual sensitivity [8x1]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
if strcmp(prob.termMode, 'rendezvous')
    dg = [de(1:6); de(15); de(8)];
else
    dg = de(1:8);
end
end

function dS = ifs_dS(Y, c)
% IFS_DS  Analytic gradient dS/dY at a node, where the switching function
% is S = 1 - ||lamV||*c/m - lamM.
%
% INPUTS:
%   Y - node augmented state [16x1]
%   c - exhaust-velocity-like constant [scalar]
% OUTPUTS:
%   dS - gradient of S w.r.t. Y [1x16]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
lamV = Y(12:14);  m = Y(7);  nv = sqrt(sum(lamV.^2));
dS = zeros(1, 16);
dS(7)      =  nv*c/m^2;                 % d/dm
dS(12:14)  = -(lamV.'/max(nv,1e-30))*c/m;
dS(15)     = -1;
end
