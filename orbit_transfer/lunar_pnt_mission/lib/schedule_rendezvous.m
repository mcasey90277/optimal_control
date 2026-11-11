function S = schedule_rendezvous(tf, sD, sA, tauD, tauT, droPhase0, slotPhase0, opts)
%% Purpose:
%
%   Assign K spacecraft waiting in a DRO to K tulip slots so the LAST
%   arrival is as early as possible, using catalog transfer cells
%   (departure phase sD(iD), arrival phase sA(iA), min-time tf(iD, iA))
%   and whole-revolution waits. A pure function of tables.
%
%   CLOCK AND PHASE CONVENTIONS (Global Constraints). Clock time t = 0
%   ND. Spacecraft k is at DRO phase mod(droPhase0(k) + t/tauD, 1); slot j
%   is at tulip phase mod(slotPhase0(j) + t/tauT, 1). Candidate for
%   spacecraft k, slot j, cell (iD, iA), wait n in 0:nMaxRevs:
%       tDep  = tauD*mod(sD(iD) - droPhase0(k), 1) + n*tauD
%       tArr  = tDep + tf(iD, iA)
%       resid = wrap(sA(iA) - mod(slotPhase0(j) + tArr/tauT, 1))
%   with wrap into (-0.5, 0.5]; admissible iff abs(resid) <= opts.tol.
%   The departure-phase difference is snapped to 0 when it lies within
%   1e-9 below 1 (floating noise in a phase equal to the DRO phase would
%   otherwise cost a whole revolution).
%
%   Assignment: minimise max(tArr) over permutations, each pair
%   represented by its earliest admissible candidate. Ties within 1e-12 go
%   to the smaller total propellant (propKgPerNd * tf of the chosen
%   cells), then to the lexicographically smaller permutation. The
%   propellant tie-break compares assignments using each pair's EARLIEST
%   candidate only: a later-but-cheaper candidate that still meets the
%   bottleneck is not considered, because the objective is the earliest
%   complete constellation and propellant is reported, not optimised.
%   Within a pair, candidates are sorted by [tArr, tf] with tArr compared
%   to 1e-12 (on a tie the earliest candidate is the cheaper one); any
%   remaining tie keeps (iD, iA, n) order.
%   Throws schedule_rendezvous:noAssignment (message lists the pairs with
%   no candidate and the K x K best-arrival table) when no permutation is
%   complete.
%
%% Inputs:
%
%  tf          nD x nA     Catalog min-time table [ND]
%  sD          1 x nD      Catalog departure phases [fraction]
%  sA          1 x nA      Catalog arrival phases [fraction]
%  tauD        scalar      DRO period [ND]
%  tauT        scalar      Tulip period [ND]
%  droPhase0   1 x K       DRO phase of spacecraft k at t = 0 [fraction]
%  slotPhase0  1 x K       Tulip phase of slot j at t = 0 [fraction]
%  opts        struct, optional
%    .tol          scalar  Arrival phase tolerance [0.005 of a tulip period]
%    .nMaxRevs     scalar  Largest wait, whole DRO revolutions [10]
%    .propKgPerNd  scalar  Propellant per ND of burn [kg/ND] [1]
%    .exclude      m x 5   Rows [k j iD iA nRev]: candidates made
%                          INADMISSIBLE (refused polishes; the study's
%                          schedule -> polish -> exclude loop) [zeros(0,5)];
%                          empty = no exclusion, bit for bit. Excluded
%                          candidates are absent from ranked and
%                          nAdmissible. Not m x 5 numeric:
%                          schedule_rendezvous:exclude
%
%% Outputs:
%
%  S.assign       1 x K    Slot of each spacecraft
%  S.iD S.iA S.nRev S.tDep S.tArr S.resid S.tf S.propKg
%                 1 x K    The chosen candidate per spacecraft
%  S.tLast        scalar   Last arrival [ND]
%  S.totalPropKg  scalar   Sum of S.propKg [kg]
%  S.best         struct   K x K matrices tArr iD iA nRev resid tDep of the
%                          earliest admissible candidate (NaN if none)
%  S.ranked       {K x K}  Each m x 6 [tArr iD iA nRev resid tDep], ALL
%                          admissible candidates sorted by tArr (polish
%                          fallback list)
%  S.nAdmissible  K x K    Admissible-candidate counts
%  S.opts         struct   Options as used
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  opts.exclude (refused-polish exclusion list)     (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    setup_lunar_pnt_paths();
    M = mission_inputs();
    o = struct('tol', 0.005, 'nMaxRevs', 10, 'propKgPerNd', M.propKgPerNd);
    t0 = tic;
    S = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, ...
                            M.droPhase0, M.slotPhase0, o);
    wall = toc(t0);
    d = S.best.tArr*M.tStar/86400;
    fprintf('Best arrival [days], rows spacecraft, columns slot:\n');
    fprintf('%s\n', sprintf([repmat(' %8.3f', 1, size(d, 2)) '\n'], d.'));
    fprintf('Admissible candidates per pair:\n');
    fprintf('%s\n', sprintf([repmat(' %5d', 1, size(d, 2)) '\n'], S.nAdmissible.'));
    fprintf('%-4s %-5s %-4s %-4s %-3s %10s %10s %10s %9s\n', 'sc', 'slot', ...
            'iD', 'iA', 'n', 'tDep[d]', 'tArr[d]', 'resid', 'prop[kg]');
    for kk = 1:numel(S.assign)
        fprintf('%-4d %-5d %-4d %-4d %-3d %10.3f %10.3f %10.5f %9.3f\n', kk, ...
                S.assign(kk), S.iD(kk), S.iA(kk), S.nRev(kk), ...
                S.tDep(kk)*M.tStar/86400, S.tArr(kk)*M.tStar/86400, ...
                S.resid(kk), S.propKg(kk));
    end
    fprintf('Last arrival %.3f d  total propellant %.3f kg  wall %.3f s\n', ...
            S.tLast*M.tStar/86400, S.totalPropKg, wall);
    return
end

if nargin < 8, opts = struct(); end
if ~isfield(opts, 'tol'),         opts.tol = 0.005;     end
if ~isfield(opts, 'nMaxRevs'),    opts.nMaxRevs = 10;   end
if ~isfield(opts, 'propKgPerNd'), opts.propKgPerNd = 1; end
if ~isfield(opts, 'exclude'),     opts.exclude = zeros(0, 5); end
ex = opts.exclude;
if isempty(ex), ex = zeros(0, 5); end
if ~isnumeric(ex) || size(ex, 2) ~= 5
    error('schedule_rendezvous:exclude', ...
          'opts.exclude must be m x 5 rows [k j iD iA nRev] (got %s)', mat2str(size(ex)));
end

K  = numel(droPhase0);
assert(numel(slotPhase0) == K, 'schedule_rendezvous:size', 'K mismatch');
sD = sD(:);  sA = sA(:).';
[nD, nA] = size(tf);
assert(numel(sD) == nD && numel(sA) == nA, 'schedule_rendezvous:size', 'table size');
W  = opts.nMaxRevs + 1;
nn = reshape(0:opts.nMaxRevs, 1, 1, W);

bestT = nan(K, K);  bestD = nan(K, K);  bestA = nan(K, K);
bestN = nan(K, K);  bestR = nan(K, K);  bestDep = nan(K, K);
bestF = nan(K, K);
nAdm  = zeros(K, K);
ranked = cell(K, K);
for kk = 1:K
    dph = mod(sD - droPhase0(kk), 1);
    dph(dph > 1 - 1e-9) = 0;
    tDep = tauD*dph + tauD*nn;                         % nD x 1 x W
    tArr = tDep + tf;                                  % nD x nA x W
    for jj = 1:K
        ph  = mod(slotPhase0(jj) + tArr/tauT, 1);
        r   = mod(sA - ph, 1);                         % [0, 1)
        r(r > 0.5) = r(r > 0.5) - 1;                   % (-0.5, 0.5]
        adm = find(abs(r) <= opts.tol);                % linear, iD fastest
        exk = ex(ex(:, 1) == kk & ex(:, 2) == jj, 3:5);
        if ~isempty(exk) && ~isempty(adm)
            [id, ia, iw] = ind2sub([nD nA W], adm);
            adm = adm(~ismember([id, ia, iw - 1], exk, 'rows'));
        end
        nAdm(kk, jj) = numel(adm);
        if isempty(adm)
            ranked{kk, jj} = zeros(0, 6);
            continue
        end
        [id, ia, iw] = ind2sub([nD nA W], adm);
        cand = [tArr(adm), id, ia, iw - 1, r(adm), tDep(sub2ind([nD 1 W], id, ones(size(id)), iw))];
        [~, ord] = sortrows([round(cand(:, 1)*1e12), tf(sub2ind([nD nA], id, ia))]);
        cand = cand(ord, :);
        ranked{kk, jj} = cand;
        bestT(kk, jj) = cand(1, 1);  bestD(kk, jj) = cand(1, 2);
        bestA(kk, jj) = cand(1, 3);  bestN(kk, jj) = cand(1, 4);
        bestR(kk, jj) = cand(1, 5);  bestDep(kk, jj) = cand(1, 6);
        bestF(kk, jj) = tf(cand(1, 2), cand(1, 3));
    end
end

S.best = struct('tArr', bestT, 'iD', bestD, 'iA', bestA, 'nRev', bestN, ...
                'resid', bestR, 'tDep', bestDep);
S.ranked = ranked;
S.nAdmissible = nAdm;
S.opts = opts;

P = sortrows(perms(1:K));                              % lexicographic
tPerm = inf(size(P, 1), 1);  prop = inf(size(P, 1), 1);
tIn = bestT;  tIn(isnan(tIn)) = inf;
pIn = opts.propKgPerNd*bestF;  pIn(isnan(pIn)) = inf;
lin = (1:K);
for kp = 1:size(P, 1)
    idx = sub2ind([K K], lin, P(kp, :));
    tPerm(kp) = max(tIn(idx));
    prop(kp)  = sum(pIn(idx));
end
if ~any(isfinite(tPerm))
    [rr, cc] = find(nAdm == 0);
    pairs = sprintf(' (sc %d, slot %d)', [rr cc].');
    tab = sprintf([repmat(' %9.4f', 1, K) ' |'], bestT.');
    error('schedule_rendezvous:noAssignment', ...
          ['no complete assignment; pairs with no admissible candidate:%s;' ...
           ' best-arrival table by row [ND]:%s'], pairs, tab);
end
tmin = min(tPerm);
c1 = find(tPerm <= tmin + 1e-12);
[~, c2] = min(prop(c1));
pmin = prop(c1(c2));
c3 = c1(prop(c1) <= pmin + 1e-12);
kp = c3(1);                                            % P is lexicographic

p = P(kp, :);
idx = sub2ind([K K], lin, p);
S.assign = p;
S.iD = bestD(idx);  S.iA = bestA(idx);  S.nRev = bestN(idx);
S.tDep = bestDep(idx);  S.tArr = bestT(idx);  S.resid = bestR(idx);
S.tf = bestF(idx);
S.propKg = opts.propKgPerNd*S.tf;
S.tLast = max(S.tArr);
S.totalPropKg = sum(S.propKg);
end
