function Q = sheet_to_catalog_file(S, ribs, outMat, opts)
%% Purpose:
%
%   Convert an arrival-phase sheet (build_arrival_sheet) plus its departure
%   ribs (rib_from_crossing) into the SHEET FILE layout that
%   build_costate_catalog_family consumes, and save it.
%
%   No schema work is needed to ship these entries: the min-time catalog
%   schema already keys every sheet by sD_frac x sA_frac with a thrust
%   axis, so a one-rung 70 mN phase sheet is a LAYOUT conversion. It does
%   need its own catalog, because the thruster differs from the shipped
%   DRO -> tulip catalog (Isp 900 s here, 1710 s there) and a catalog
%   carries one thruster.
%
%   Only CERTIFIED entries are placed: OK(iD,iA) is false and t_f is NaN
%   wherever nothing passed the gate stack, which is what makes the
%   catalog's availability grid honest.
%
%% Inputs:
%
%  S                        struct                  arrival sheet: .sA
%                                                   [1 x nA], .TF (days),
%                                                   .Z8 [8 x nA], .cand
%  ribs                     cell | []               rib_from_crossing
%                                                   outputs (.pts with .sD
%                                                   .sA .z .tfDays)
%  outMat                   char                    sheet-file path
%  opts                     struct (optional)
%   .nD [12] .sD0 [0] .thrustN [0.070] .ispS [900] .m0kg [150]
%   .tauDRO [1] .NpTulip [7] .pmTulip [-1]
%
%% Outputs:
%
%  Q                        struct                  .OK [nD x nA x 1]
%                                                   .Z8 [8 x nD x nA x 1]
%                                                   .TF [nD x nA x 1] ND
%                                                   .sD .sA .rungs .meta
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
nD = d('nD', 12);  sD0 = d('sD0', 0);  nA = numel(S.sA);
tStar = 382981.289129055;  lStar = 389703.264829278;

Q = struct();
Q.sD = mod(sD0 + (0:nD-1)/nD, 1);
Q.sA = S.sA(:)';
Q.rungs = d('thrustN', 0.070);
Q.OK = false(nD, nA, 1);
Q.TF = nan(nD, nA, 1);
Q.Z8 = nan(8, nD, nA, 1);

% ---- the spine: the certified minimum at each arrival phase, at sD0 -----
iD0 = idxOf(Q.sD, sD0);
for j = 1:nA
    if ~isfinite(S.TF(j)), continue, end
    Q.OK(iD0, j, 1) = true;
    Q.TF(iD0, j, 1) = S.TF(j)*86400/tStar;          % days -> ND
    Q.Z8(:, iD0, j, 1) = S.Z8(:, j);
end

% ---- the ribs: certified departure points off the spine ----------------
if nargin >= 2 && ~isempty(ribs)
    if ~iscell(ribs), ribs = {ribs}; end
    for k = 1:numel(ribs)
        R = ribs{k};
        if ~isfield(R, 'pts') || isempty(R.pts), continue, end
        for m = 1:numel(R.pts)
            P = R.pts(m);
            if ~P.ok, continue, end
            iD = idxOf(Q.sD, P.sD);  iA = idxOf(Q.sA, P.sA);
            assert(~isempty(iD) && ~isempty(iA), 'rib point (%.4f, %.4f) is off the grid', P.sD, P.sA);
            if Q.OK(iD, iA, 1) && Q.TF(iD, iA, 1) <= P.z(8), continue, end   % keep the faster
            Q.OK(iD, iA, 1) = true;
            Q.TF(iD, iA, 1) = P.z(8);
            Q.Z8(:, iD, iA, 1) = P.z(:);
        end
    end
end

% ---- meta: what the packager reads -------------------------------------
tauDRO = d('tauDRO', 1.0);  Np = d('NpTulip', 7);  pm = d('pmTulip', -1);
Q.meta = struct('muStar', 0.012150585609624, 'lStar', lStar, 'tStar', tStar, ...
    'ispS', d('ispS', 900), 'm0kg', d('m0kg', 150), 'tauDRO', tauDRO, ...
    'depFamily', 'dro', 'depParams', struct('tau', tauDRO), ...
    'arrFamily', 'tulip', 'arrParams', struct('Np', Np, 'pm', pm, 'tau', 5*2*pi/6), ...
    'NpTulip', Np, 'pmTulip', pm, 'periodTulip', 5*2*pi/6);

if ~isempty(outMat), save(outMat, '-struct', 'Q'); end
end

function k = idxOf(grid, v)
% IDXOF  Index of the grid entry matching phase v (mod 1).  INPUTS: grid; v.
% OUTPUTS: k (empty if none).
k = find(abs(mod(grid - v + 0.5, 1) - 0.5) < 1e-8, 1);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
