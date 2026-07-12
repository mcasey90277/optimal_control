function recon = ifs_reconstruct(Z, prob, nptsTotal)
% IFS_RECONSTRUCT  Integrate an IFS solution into a seed-layout trajectory.
%
% An IFS unknown vector Z is initial costates + interior node states + switch
% times -- not a sampled trajectory. This propagates the k+1 arcs of the
% 16-dim PMP system (IFS_EOM, hard throttle uArc per arc) and assembles the
% state, control, and costate histories on a fine mesh, in the SAME layout the
% PSR data products / movie use (out.X [8xnN] = [r;v;m;t], out.U [4xnN] =
% [alpha;s]). This lets an IFS result be exported and animated with the same
% tooling as a direct/PSR solution.
%
% The thrust DIRECTION is the primer alpha = -lamV/||lamV|| on burn arcs (the
% actual PMP control, recovered from the costates -- unlike the direct solver,
% where it is a free NLP variable) and is left as the propagated value on coast
% arcs (thrust off, direction irrelevant).
%
% INPUTS:
%   Z         - IFS unknown vector [(8+17k)x1]
%   prob      - IFS problem struct (rv0, m0, t0, tau0, tauf, k, uArc, Tmax, c,
%               muStar, pSund, odeOpts, tauParam, rvf, factor if present)
%   nptsTotal - target number of mesh points across all arcs [default 4000]
%
% OUTPUTS:
%   recon - struct in seed layout + extras:
%     .out.X [8xnN], .out.U [4xnN], .out.mf, .out.switches, .out.maxDefect(NaN)
%     .sigma [nNx1] normalized (tau/tauf), .tauf0 (=prob.tauf), .tau [1xnN]
%     .rv0 [1x6], .rvf [1x6], .factor
%     .lam [8xnN] costates [lamR;lamV;lamM;lamT]
%     .S   [1xnN] switching function 1 - ||lamV||c/m - lamM
%     .tauSwitch [1xk]
%
% REFERENCES: ifs_certify.m (arc-propagation pattern), PSR/psr_movie.m (layout).

if nargin < 3 || isempty(nptsTotal), nptsTotal = 4000; end
k = prob.k;
if ~isfield(prob,'tauParam'), prob.tauParam = 'sigmoid'; end

[lam0, N, gblk] = ifs_unpack(Z, k);
tau   = ifs_taus(gblk, prob.tau0, prob.tauf, prob.tauParam);
edges = [prob.tau0, tau(:).', prob.tauf];               % k+2 arc boundaries

% arc-start augmented states: arc 1 from (rv0; costate), arcs 2..k+1 from nodes
startY = cell(1, k+1);
startY{1} = [prob.rv0(:); prob.m0; prob.t0; lam0];
for a = 2:k+1, startY{a} = N(:, a-1); end

% distribute mesh points across arcs by arc length (>= 8 per arc)
arcLen = diff(edges);
npa = max(8, round(nptsTotal * arcLen / max(sum(arcLen), eps)));

Y = [];  tauAll = [];
for a = 1:k+1
    sp = [edges(a), edges(a+1)];
    if sp(2) - sp(1) <= 1e-13, continue; end            % collapsed gap: skip
    tgrid = linspace(sp(1), sp(2), npa(a));
    [~, Ya] = ode113(@(s,y) ifs_eom(s, y, prob.Tmax, prob.c, prob.muStar, ...
                     prob.pSund, prob.uArc(a)), tgrid, startY{a}, prob.odeOpts);
    if a > 1, Ya = Ya(2:end, :);  tgrid = tgrid(2:end);  end   % drop duplicate joint
    Y = [Y; Ya];  tauAll = [tauAll, tgrid]; %#ok<AGROW>
end
Y = Y.';                                                % 16 x nN
nN = size(Y, 2);

X   = Y(1:8, :);                                        % [r;v;m;t]
lam = Y(9:16, :);                                       % [lamR;lamV;lamM;lamT]
lamV = lam(4:6, :);
nrm  = sqrt(sum(lamV.^2, 1));
alpha = -lamV ./ max(nrm, 1e-12);                       % primer direction
% throttle per node = its arc's uArc (map by which arc each tau falls in)
sthr = zeros(1, nN);
for j = 1:nN
    a = find(tauAll(j) >= edges(1:end-1) & tauAll(j) <= edges(2:end), 1, 'last');
    if isempty(a), a = k+1; end
    sthr(j) = prob.uArc(a);
end
U = [alpha; sthr];
S = 1 - nrm.*prob.c./X(7,:) - lam(7,:);

recon.out = struct('X', X, 'U', U, 'mf', X(7,end), ...
    'switches', sum(abs(diff(sthr > 0.5))), 'maxDefect', NaN);
recon.sigma = (tauAll(:) / prob.tauf);                 % normalized [0,1]
recon.tauf0 = prob.tauf;
recon.tau   = tauAll;
recon.rv0   = prob.rv0(:).';
if isfield(prob,'rvf'), recon.rvf = prob.rvf(:).'; else, recon.rvf = X(1:6,end).'; end
if isfield(prob,'factor'), recon.factor = prob.factor; else, recon.factor = NaN; end
recon.lam = lam;  recon.S = S;  recon.tauSwitch = tau(:).';
end
