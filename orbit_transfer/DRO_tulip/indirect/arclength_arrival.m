function varargout = arclength_arrival(arg, opts)
%% Purpose:
%
%   ARRIVAL-PHASE pseudo-arclength continuation of the minimum-time
%   DRO -> tulip transfer at fixed thrust and departure phase, in the
%   homogeneous chart (rho free on the sphere), on the generic engine
%   arclength_ms. This is the sheet's spine: it follows the solution curve
%   through folds and records every crossing of the requested arrival-phase
%   grid levels -- repeated crossings included, so a grid point can hold
%   more than one candidate. Certification and choice happen downstream.
%
%   Two calls:
%     [B, anc] = arclength_arrival('setup', opts)   problem closures + anchor
%     A        = arclength_arrival(anc, opts)       trace an arc from anc
%
%   Why arrival phase is a cheap parameter here (GPT-6 Astra, 2026-09-09):
%   it enters ONLY the six terminal state-matching rows,
%       R_sA = [0; -x_A'(sA); 0; 0; 0],   x_A'(sA) = T_A * f_cr3bp(x_A(sA)),
%   so the derivative is analytic (cr3bp_field) -- no finite differencing
%   of the whole residual, unlike thrust.
%
%   Scaling is per BLOCK, as the engine demands: states 1, costates 1 (they
%   live on the unit sphere in this chart), t_f by its anchor value, rho 1,
%   phase 1.
%
%% Inputs:
%
%  arg                      'setup' | anchor struct
%
%  opts                     struct (optional)
%   setup:  .thrustN [0.070] .ispS [900] .m0kg [150] .tauDRO [1] .NpTulip [7]
%           .sD [0] departure phase, .anchorMat [results/mintime_70mN_anchor.mat]
%           (root as z + it.Y, or best.z + best.it.Y), .sA0 [0.0754] the
%           anchor's arrival phase, .K [from the anchor]
%   arc:    .direction [+1] .sAStop [inf] (stop beyond) .ds [0.002]
%           .dsMin [1e-5] .dsMax [0.02] .nStep [400] .levels []
%           .deadlineSec [inf] .logFile ''
%
%% Outputs:
%
%  setup: B (closures: .stateA .stateD .tauA .tauD .mu .Tnd .cnd .res(sA)
%            .dRdq(p,sA) .Dx .rv0), anc (.p .sA .sD .ctf .nExtra .termRows .K)
%  arc:   A (arclength_ms output; A.p{k} are full ms unknown vectors)
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));

if ischar(arg) && strcmp(arg, 'setup')
    thrustN = d('thrustN', 0.070);  ispS = d('ispS', 900);  m0kg = d('m0kg', 150);
    ob = struct('muStar', 0.012150585609624, 'lStar', 389703.264829278, ...
                'tStar', 382981.289129055, 'tauDRO', d('tauDRO', 1.0), ...
                'NpTulip', d('NpTulip', 7), 'tauTulip', 5*2*pi/6, 'pmTulip', -1, ...
                'ispS', ispS, 'm0kg', m0kg);
    lStar = ob.lStar;  tStar = ob.tStar;  mu = ob.muStar;
    g0  = 9.80665*tStar^2/(1000*lStar);
    cnd = (ispS/tStar)*g0;
    Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);
    [tD, rvD, tT, rvT] = ladder_endpoints(ob);
    B = struct();
    B.mu = mu;  B.Tnd = Tnd;  B.cnd = cnd;  B.tauA = tT(end);  B.tauD = tD(end);
    B.stateD = @(s) interp1(tD, rvD, mod(s,1)*tD(end), 'spline')';
    B.stateA = @(s) interp1(tT, rvT, mod(s,1)*tT(end), 'spline')';
    sD = d('sD', 0);
    B.rv0 = B.stateD(sD);

    % the anchor: a certified rho = 1 root, re-normalised onto the sphere
    anchorMat = d('anchorMat', fullfile(here, 'results', 'mintime_70mN_anchor.mat'));
    Aanc = load(anchorMat);
    % two storage layouts: the demo anchor (z, it.Y, Tnd, cnd) and the
    % certified sheet cells (best.z, best.it.Y, no engine constants)
    if isfield(Aanc, 'best'), root = Aanc.best; else, root = Aanc; end
    z = root.z(:);  Y = root.it.Y;  K = d('K', size(Y, 2));
    if K ~= size(Y, 2)
        % A DIFFERENT MESH than the stored root's: re-cut the same flight
        % into K junctions (seed_from_z8) rather than pad or drop columns.
        % Used by the mesh-independence check on the arclength metric.
        sd0 = seed_from_z8(z, B.rv0(1:6), K, Tnd, cnd, mu);
        Y = sd0.Y(:, 1:K);
    end
    if isfield(Aanc, 'Tnd') && isfield(Aanc, 'cnd')
        assert(abs(Aanc.Tnd - Tnd)/Tnd < 1e-10 && abs(Aanc.cnd - cnd)/cnd < 1e-10, ...
               'anchor certified at a different operating point');
    end
    sA0 = d('sA0', 0.0754);          % the anchor's grid phase (caller's claim)
    seed = struct('tf', z(8), 'tGrid', linspace(0, z(8), K+1), 'Y', [Y, Y(:,end)]);
    seed.Y(1:7, 1) = [B.rv0(1:6); 1];
    % RE-SOLVE the anchor at THIS grid's exact arrival state rather than
    % trusting the stored root: the certificate was computed against an
    % arrival state 49 km from stateA(0.0754) (a rounded phase), and a
    % stored root is only a root of the problem it was solved for. One
    % Newton polish absorbs that; the assert below is the real check.
    [zh, ih] = ms_tfmin_hom(B.rv0(1:6), B.stateA(sA0), seed, Tnd, cnd, mu, ...
                            struct('tolR', 1e-11, 'wallSec', 120));
    assert(ih.converged, 'anchor did not re-converge at the grid phase (|R| = %.1e)', ih.normR);
    % the re-solve may only absorb phase rounding (~1e-4 in t_f). A larger
    % move means the stored root belongs to another operating point or
    % phase -- the guard that replaces the Tnd/cnd assert when the .mat
    % carries no engine constants.
    assert(abs(zh(8) - z(8))/z(8) < 1e-3, ...
           'anchor moved %.2e in t_f on re-solve: wrong phase or operating point', abs(zh(8) - z(8))/z(8));
    seed.Y = [ih.Y, ih.Y(:,end)];  seed.tf = zh(8);  seed.extra = ih.rho;
    seed.tGrid = linspace(0, zh(8), K+1);
    [~, info] = ms_tfmin_hom(B.rv0(1:6), B.stateA(sA0), seed, Tnd, cnd, mu, ...
                             struct('assembleOnly', true));
    assert(norm(info.R, inf) < 1e-8, 'anchor is not a root after re-solve (|R| = %.1e)', norm(info.R, inf));
    n = numel(info.p);  ctf = info.ctf;  nX = info.nExtra;
    termRows = (K-1)*14 + (1:8);

    % residual factory in sA (the arrival state is the only thing that moves)
    B.res  = @(sA) resHandle(B.rv0, B.stateA(sA), seed, Tnd, cnd, mu);
    B.dRdq = @(p, sA) dRdsA(sA, B, n, termRows);
    % BLOCK SCALES, MESH-INDEPENDENT. The trajectory block is K junctions x
    % 14 components; at unit weight its 336 coordinates swamp the single
    % phase coordinate and an arclength step of 0.002 moves sA by 5e-5
    % (measured: 60 steps advanced 0.0015 of a period). Weighting every
    % junction coordinate by sqrt(K) makes the block's contribution
    % quadrature-like -- one junction's worth, independent of K (Astra,
    % 2026-09-09) -- so phase and trajectory motion are commensurate.
    Dx = sqrt(K)*ones(n, 1);
    Dx(ctf) = zh(8);                 % t_f by its own value
    Dx(end) = 1;                     % rho: O(0.01-0.1), keep at unit scale
    B.Dx = Dx;
    anc = struct('p', info.p, 'sA', sA0, 'sD', sD, 'ctf', ctf, 'nExtra', nX, ...
                 'termRows', termRows, 'K', K, 'n', n);
    varargout = {B, anc};
    return
end

% ---- arc mode --------------------------------------------------------
anc = arg;
[B, ~] = arclength_arrival('setup', opts);
lvl = d('levels', []);
sAStop = d('sAStop', inf);
dirn = d('direction', +1);
if dirn > 0, qStop = [-inf sAStop]; else, qStop = [sAStop inf]; end
adm = @(p, q) p(anc.ctf) > 0 && all(isfinite(p));
A = arclength_ms(B.res, B.dRdq, anc.p, anc.sA, struct( ...
    'direction', dirn, 'Dx', B.Dx, 'sq', 1, ...
    'ds', d('ds', 0.01), 'dsMin', d('dsMin', 1e-5), 'dsMax', d('dsMax', 0.2), ...
    'nStep', d('nStep', 400), 'qStop', qStop, 'levels', lvl, ...
    'newtonTol', d('newtonTol', 1e-9), 'admissible', adm, ...
    'deadlineSec', d('deadlineSec', inf), 'logFile', d('logFile', '')));
A.anc = anc;  A.B = rmfield(B, {'res', 'dRdq', 'stateA', 'stateD'});
varargout = {A};
end

% ------------------------------------------------------------------------
function h = resHandle(rv0, rvf, seed, Tnd, cnd, mu)
% RESHANDLE  Production homogeneous ms residual for one arrival state.
% INPUTS: rv0; rvf; seed; Tnd; cnd; mu.  OUTPUTS: h (p -> [R, J]).
sd = seed;  sd.Y(1:7, 1) = [rv0(1:6); 1];
[~, inf_] = ms_tfmin_hom(rv0(1:6), rvf(1:6), sd, Tnd, cnd, mu, ...
    struct('assembleOnly', true, 'handleOnly', true, 'tfLo', 0.05, 'tfHi', 20, 'pMax', 1e7));
h = inf_.residual;
end

function col = dRdsA(sA, B, n, termRows)
% DRDSA  Analytic d R / d sA: the six terminal state rows carry
% -x_A'(sA) = -T_A f(x_A(sA)); everything else is zero.
% INPUTS: sA; B; n; termRows.  OUTPUTS: col [n x 1].
col = zeros(n, 1);
xA = B.stateA(sA);
col(termRows(1:6)) = -B.tauA * cr3bp_field(xA, B.mu);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
