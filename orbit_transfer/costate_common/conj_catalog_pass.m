function S = conj_catalog_pass(catMat, opts)
% CONJ_CATALOG_PASS  Run the conjugate-point test over every entry of a
%   compact costate catalog and record the verdicts.
%
%   For each stored entry: rebuild the endpoint states from the sheet's
%   family recipes (same get_family_orbit + spline construction the ladder
%   engines used), fly the stored z8 to build a K-junction seed, re-solve
%   with ms_tfmin(conjTest) seeded AT the converged solution (1 Newton
%   iteration expected), and record the free-time quotiented Jacobi verdict
%   from ms_conjugate_test. A verdict is only recorded when the re-solve
%   reproduces the stored solution (converged, |z - z8| < tolDz); otherwise
%   the entry is marked NOT RUN (-1), never FAIL.
%
%   Campaign contract (matlab-campaign discipline): file logging, sidecar
%   progress .mat saved after EVERY entry, attempt counter written BEFORE
%   each solve, clean batch-budget exit between entries, resume for free.
%   Verdicts are written back into the catalog .mat only by an explicit
%   'writeback' call after a complete census (a .bak copy is made first).
%
% INPUTS:
%   catMat - path to a catalog .mat (single variable, schema v1/v2) [char]
%   opts   - (optional) struct:
%            .logFile    [''] append-mode log (empty = stdout)
%            .batchSec   [inf] clean-exit wall budget for this call
%            .maxEntries [inf] process at most N entries this call
%            .K          [24] multiple-shooting segment count
%            .maxAtt     [2] attempts per entry before it is retired
%            .tolDz      [1e-6] |z - z8| gate for a valid re-solve
%            .wallSec    [60] per-entry ms_tfmin wall budget (advisory)
%            .sideMat    [<catMat minus .mat>_conjprog.mat] sidecar path
%            .writeback  [false] after a COMPLETE census, write conj_pass /
%                        conj_ncross / conj_atfinal grids + .conj_test meta
%                        into the catalog .mat (backs up to .bak_conj first)
%
% OUTPUTS:
%   S      - struct: .done (all entries decided), .nPass, .nFail, .nNotrun,
%            .nTodo (still undecided and not retired), .sideMat
%
% REFERENCES:
%   [1] costate_common/ms_conjugate_test.m (the instrument; header = math)
%   [2] costate_common/golden_cells.m (the per-entry reconstruction pattern)
%   [3] orbit_transfer/STATUS_AND_ROADMAP.md §6 step 1 (why this exists)

if nargin < 2, opts = struct(); end
K        = fieldd(opts, 'K', 24);
batchSec = fieldd(opts, 'batchSec', inf);
maxEnt   = fieldd(opts, 'maxEntries', inf);
maxAtt   = fieldd(opts, 'maxAtt', 2);
tolDz    = fieldd(opts, 'tolDz', 1e-6);
wallSec  = fieldd(opts, 'wallSec', 60);
logFile  = fieldd(opts, 'logFile', '');
[pth, base] = fileparts(catMat);
sideMat  = fieldd(opts, 'sideMat', fullfile(pth, [base '_conjprog.mat']));
writeback = fieldd(opts, 'writeback', false);

lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

%% Load catalog (name-agnostic single variable):
L = load(catMat);
fn = fieldnames(L);
assert(numel(fn) == 1, 'expected one variable in %s', catMat);
cat_ = L.(fn{1});
nS = numel(cat_.sheets);

%% Sidecar: per-sheet verdict + attempt grids (resume for free):
if exist(sideMat, 'file')
    P = load(sideMat);
    assert(P.K == K && numel(P.CONJ) == nS, ...
        'sidecar %s incompatible (K or sheet count changed)', sideMat);
else
    P = struct('K', K, 'catMat', catMat, 'created', datestr(now), ...
               'CONJ', {cell(1,nS)}, 'NCROSS', {cell(1,nS)}, ...
               'ATFIN', {cell(1,nS)}, 'ATT', {cell(1,nS)}, ...
               'DZ', {cell(1,nS)});
    for ks = 1:nS
        sz = size(cat_.sheets(ks).has_solution);
        P.CONJ{ks}   = -ones(sz, 'int8');   % 1 pass / 0 fail / -1 not run
        P.NCROSS{ks} = -ones(sz, 'int8');
        P.ATFIN{ks}  = -ones(sz, 'int8');
        P.ATT{ks}    = zeros(sz, 'uint8');
        P.DZ{ks}     = nan(sz);
    end
    save(sideMat, '-struct', 'P');
end

%% Physics constants (identical derivation to the deliverable examples):
mu    = cat_.constants.muStar;
lStar = cat_.constants.lStar_km;
tStar = cat_.constants.tStar_s;
cnd   = cat_.thruster.c_nd;
m0    = cat_.thruster.m0_kg;
ndT   = @(TN) (TN/m0)*tStar^2/(lStar*1000);

tAll = tic;
nDone = 0;
for ks = 1:nS
    sh = cat_.sheets(ks);
    OKg = sh.has_solution;
    todo = OKg & (P.CONJ{ks} == -1) & (P.ATT{ks} < maxAtt);
    if ~any(todo(:)), continue, end

    % Rebuild both orbits ONCE per sheet. v2 sheets carry full recipes;
    % v1 sheets (DRO catalog: none; DPO catalog: departure only) fall back
    % to the original convention: dep = dro(tauDRO), arr = tulip(Np, pm).
    if isfield(sh, 'dep_family')
        depFam = sh.dep_family;  depPar = sh.dep_params;
    else
        depFam = 'dro';          depPar = struct('tau', sh.tauDRO);
    end
    if isfield(sh, 'arr_family')
        arrFam = sh.arr_family;  arrPar = sh.arr_params;
    else
        arrFam = 'tulip';        arrPar = struct('Np', sh.Np, 'pm', sh.pm);
    end
    [tD, rvD] = get_family_orbit(depFam, depPar);
    [tA, rvA] = get_family_orbit(arrFam, arrPar);
    lg('[sheet %d/%d] %s(%g) -> %s: %d entries to test', ks, nS, ...
       depFam, sh.tauDRO, arrFam, nnz(todo));

    [nD, nA, nR] = size(OKg);
    for iD = 1:nD
     for iA = 1:nA
      for kr = 1:nR
        if ~todo(iD,iA,kr), continue, end
        if toc(tAll) > batchSec || nDone >= maxEnt
            save(sideMat, '-struct', 'P');
            lg('[batch] clean exit: %d entries this call', nDone);
            S = census(P, cat_, sideMat, maxAtt); return
        end
        z8 = sh.z8(:, sh.entry_index(iD,iA,kr));
        % Sanity screen on DISCRETE data before any integration:
        if ~all(isfinite(z8)) || z8(8) <= 0
            P.CONJ{ks}(iD,iA,kr) = -1;  P.ATT{ks}(iD,iA,kr) = maxAtt;
            save(sideMat, '-struct', 'P');
            lg('[s%d %d,%d,r%d] SKIP: bad z8', ks, iD, iA, kr); continue
        end
        % Attempt counter BEFORE solving (hang-proof resume):
        P.ATT{ks}(iD,iA,kr) = P.ATT{ks}(iD,iA,kr) + 1;
        save(sideMat, '-struct', 'P');

        rv0 = interp1(tD, rvD, mod(sh.sD_frac(iD),1)*tD(end), 'spline');
        rvf = interp1(tA, rvA, mod(sh.sA_frac(iA),1)*tA(end), 'spline');
        Tnd = ndT(cat_.rungs_N(kr));
        try
            seed = seed_from_z8(z8, rv0(1:6), K, Tnd, cnd, mu);
            [z, info] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, mu, ...
                struct('conjTest', true, 'wallSec', wallSec));
            dz = norm(z - z8);
            P.DZ{ks}(iD,iA,kr) = dz;
            if info.converged && dz < tolDz
                P.CONJ{ks}(iD,iA,kr)   = int8(info.conj.pass);
                P.NCROSS{ks}(iD,iA,kr) = int8(min(info.conj.nCrossings, 127));
                P.ATFIN{ks}(iD,iA,kr)  = int8(info.conj.atFinal);
                if ~info.conj.pass
                    lg('[s%d %d,%d,r%d] CONJUGATE POINT: nCross=%d atFinal=%d dz=%.1e', ...
                       ks, iD, iA, kr, info.conj.nCrossings, info.conj.atFinal, dz);
                end
            else
                lg('[s%d %d,%d,r%d] UNVERIFIED: conv=%d dz=%.1e (att %d/%d)', ...
                   ks, iD, iA, kr, info.converged, dz, P.ATT{ks}(iD,iA,kr), maxAtt);
            end
        catch err
            lg('[s%d %d,%d,r%d] ERROR: %s (att %d/%d)', ks, iD, iA, kr, ...
               err.message, P.ATT{ks}(iD,iA,kr), maxAtt);
        end
        save(sideMat, '-struct', 'P');
        nDone = nDone + 1;
        if mod(nDone, 50) == 0
            lg('[progress] %d entries this call, %.1f s elapsed', nDone, toc(tAll));
        end
      end
     end
    end
end

S = census(P, cat_, sideMat, maxAtt);
lg('[census] pass %d / fail %d / notrun %d / todo %d  (done=%d)', ...
   S.nPass, S.nFail, S.nNotrun, S.nTodo, S.done);

%% Writeback (explicit, only on a complete census):
if writeback
    assert(S.done, 'census incomplete (%d todo) -- not writing back', S.nTodo);
    bak = [catMat '.bak_conj'];
    if ~exist(bak, 'file'), copyfile(catMat, bak); end
    for ks = 1:nS
        cat_.sheets(ks).conj_pass    = P.CONJ{ks};
        cat_.sheets(ks).conj_ncross  = P.NCROSS{ks};
        cat_.sheets(ks).conj_atfinal = P.ATFIN{ks};
    end
    cat_.conj_test = struct('date', datestr(now, 'yyyy-mm-dd'), ...
        'instrument', 'costate_common/ms_conjugate_test (free-time quotiented Jacobi)', ...
        'K', K, 'tolDz', tolDz, 'nPass', S.nPass, 'nFail', S.nFail, ...
        'nNotrun', S.nNotrun, 'meaning', ['conj_pass: 1 = no conjugate point in ' ...
        '(0, tf) sampled at K-1 junctions (necessary condition PASSES); ' ...
        '0 = sign change strictly inside (CONJUGATE POINT: not a local min); ' ...
        '-1 = not run / re-solve did not reproduce the entry']);
    Lout = struct(fn{1}, cat_);
    save(catMat, '-struct', 'Lout');
    lg('[writeback] verdicts stored in %s (backup: %s)', catMat, bak);
end
end

% ------------------------------------------------------------------------
function S = census(P, cat_, sideMat, maxAtt)
% CENSUS  Count verdicts from the DATA (sidecar), not the log.
nPass = 0; nFail = 0; nNotrun = 0; nTodo = 0;
for ks = 1:numel(cat_.sheets)
    OKg = cat_.sheets(ks).has_solution;
    C = P.CONJ{ks}(OKg); A = P.ATT{ks}(OKg);
    nPass = nPass + nnz(C == 1);
    nFail = nFail + nnz(C == 0);
    nNotrun = nNotrun + nnz(C == -1 & A >= maxAtt);
    nTodo = nTodo + nnz(C == -1 & A < maxAtt);
end
S = struct('done', nTodo == 0, 'nPass', nPass, 'nFail', nFail, ...
           'nNotrun', nNotrun, 'nTodo', nTodo, 'sideMat', sideMat);
end

% ------------------------------------------------------------------------
function logmsg(f, s)
% LOGMSG  Append to file (or stdout if none).  INPUTS: f;s.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s %s\n', datestr(now,'HH:MM:SS'), s); fclose(fid);
end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
