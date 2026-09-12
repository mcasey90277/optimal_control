function S = second_order_pass(catMat, opts)
%% Purpose:
%
%   Sweep the three 2026-09-10 second-order instruments over every entry of a
%   compact costate catalog and record what they measure, so the caveats
%   become NUMBERS PER ENTRY instead of prose in a document:
%
%     conj_spectrum  dense singular-SPECTRUM scan -- interior crossings and
%                    multiplicity, closing the sampled sign test's two blind
%                    spots (two crossings in one segment; an even-order zero)
%     lift_margin    the dim S rank statement as an Eckart-Young margin
%                    against a MEASURED error, not a threshold -- C is built
%                    at two integration tolerances and their difference IS
%                    the error
%     h6_margin      lambda_m(0) < c/T, excluding the reduced problem's
%                    spurious-zero mechanism
%
%   Campaign discipline, as `conj_catalog_pass`: file logging, a sidecar
%   progress .mat written after EVERY entry, a clean wall-budget exit between
%   entries, resume for free, and writeback into the catalog only on an
%   explicit call after a COMPLETE census.
%
%   Nothing here can turn a stored entry into a failure on its own: these
%   measure, and a measurement that comes out badly is reported, not acted
%   on. Deciding what a bad margin means is a separate, human step.
%
%% Inputs:
%
%  catMat                   char                    catalog .mat
%  opts                     struct (optional)
%   .logFile [''] .sideMat [<catMat>_2ndprog.mat] .batchSec [inf]
%   .maxEntries [inf] .K [24] .nSub [8] .writeback [false]
%   .relTolPair [1e-12 1e-9] the two integration settings whose difference
%   is the lift matrix's MEASURED error. The first sweep used [1e-10 1e-7]
%   and read 4-9x on seven 26-day entries whose sigma_6 is 0.99: the loose
%   setting's own error was the whole estimate (43x at this pair, 718x at
%   [1e-12 1e-10] on entry (1,10), FINDINGS 42)
%   .adoptLegacy [false] accept a sidecar written before records carried
%   their entry identity, ONLY if every written-back value in the catalog
%   equals the sidecar's; its records are then keyed and saved
%
%% Outputs:
%
%  S                        struct                  .done .nDone .nTodo
%                                                   .worstSpectrum .worstH6
%                                                   .worstLift .rows
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  09/11/2026  every sidecar record carries the IDENTITY of the entry it
%              measured (cell key + z8) and a resumed record must match the
%              catalog -- the sidecar was positional, and the chain script
%              pointed at the first sweep's sidecar, whose lift margins
%              differ from the catalog's by up to 1.35e4
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
sideMat = d('sideMat', [strrep(catMat, '.mat', '') '_2ndprog.mat']);
logFile = d('logFile', '');  batchSec = d('batchSec', inf);
maxEnt = d('maxEntries', inf);  K = d('K', 24);  nSub = d('nSub', 8);
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
t0 = tic;

L = load(catMat);  fn = fieldnames(L);  cat_ = L.(fn{1});
assert(numel(cat_.sheets) == 1, 'this pass handles a single-sheet phase catalog');
s = cat_.sheets(1);
lStar = cat_.constants.lStar_km;  tStar = cat_.constants.tStar_s;  mu = cat_.constants.muStar;
m0 = cat_.thruster.m0_kg;  isp = cat_.thruster.isp_s;
% THE shared propulsion conversion (costate_common/nd_propulsion)
ndp = nd_propulsion(cat_.rungs_N(1), isp, m0, lStar, tStar);
cnd = ndp.cnd;   Tnd = ndp.Tnd;

ob = struct('muStar', mu, 'lStar', lStar, 'tStar', tStar, 'tauDRO', s.tauDRO, ...
            'NpTulip', s.Np, 'tauTulip', s.period_tulip_nd, 'pmTulip', s.pm, ...
            'ispS', isp, 'm0kg', m0);
[tD, rvD, ~, ~] = ladder_endpoints(ob);
stD = phase_state(tD, rvD);        % THE shared endpoint rule (FINDINGS 44)

[iD, iA, iR] = ind2sub(size(s.has_solution), find(s.has_solution));
n = numel(iD);
% Every record carries the IDENTITY of the entry it measured -- its cell
% (iD, iA, iR) and its z8 -- and a resumed record must match the catalog it
% is resumed against. A positional sidecar resumed against a re-packaged
% catalog with a different entry set would hand one entry's measurements to
% another without a sound.
keys = [iD(:) iA(:) iR(:)];
Z = zeros(8, n);
for q = 1:n, Z(:, q) = s.z8(:, s.entry_index(iD(q), iA(q), iR(q))); end
if isfile(sideMat)
    P = load(sideMat);  R = P.R;
    if numel(R) ~= n
        error('second_order_pass:staleSidecar', ...
              'sidecar %s holds %d records but the catalog has %d entries', sideMat, numel(R), n);
    end
    if ~isfield(R, 'key')
        if ~d('adoptLegacy', false)
            error('second_order_pass:unkeyedSidecar', ...
                  ['sidecar %s carries no entry identity (written before 2026-09-11). Pass ' ...
                   'adoptLegacy = true to adopt it if it matches the catalog''s written-back ' ...
                   'values, or name a fresh sidecar'], sideMat);
        end
        legacyMatches(R, s, keys, sideMat);
        for q = 1:n, R(q).key = keys(q,:);  R(q).z8 = Z(:,q); end
        save(sideMat, 'R');
        lg('adopted legacy sidecar %s: every written-back value matches the catalog; records keyed', sideMat);
    end
    for q = 1:n
        if ~R(q).done
            R(q).key = keys(q,:);  R(q).z8 = Z(:,q);       % not measured yet: bind to this build
        elseif ~(isequal(R(q).key, keys(q,:)) && isequal(R(q).z8(:), Z(:,q)))
            % sprintf first: error() refuses the non-scalar key arguments
            msg = sprintf(['record %d of %s measured cell (%d,%d,%d), but position %d of this ' ...
                           'catalog is cell (%d,%d,%d)%s: the sidecar belongs to a different ' ...
                           'catalog build'], q, sideMat, R(q).key, q, keys(q,:), ...
                          tern(isequal(R(q).key, keys(q,:)), ' with a different z8', ''));
            error('second_order_pass:staleSidecar', '%s', msg);
        end
    end
    lg('resuming from %s (%d of %d already measured, identities checked)', sideMat, nnz([R.done]), n);
else
    R = struct('done', num2cell(false(1,n)), 'nInterior', [], 'multiplicity', [], ...
               'minRelSigma', [], 'nInteriorCand', [], 'nNearMiss', [], 'nZero', [], ...
               'nUnresolved', [], 'nEndCand', [], 'conjClear', [], ...
               'candidates', [], 'h6margin', [], 'h6ok', [], 'h6clearance', [], ...
               'liftMargin', [], 'liftCertified', [], 'relTolPair', [], 'key', [], 'z8', []);
    R = R(:).';
    for q = 1:n, R(q).key = keys(q,:);  R(q).z8 = Z(:,q); end
end

relPair = d('relTolPair', [1e-12 1e-9]);
nThis = 0;
for q = 1:n
    if R(q).done, continue, end
    if toc(t0) > batchSec, lg('batch budget reached, exiting cleanly'); break, end
    if nThis >= maxEnt, lg('maxEntries reached'); break, end
    z8 = s.z8(:, s.entry_index(iD(q), iA(q), iR(q)));
    rv0 = stD(s.sD_frac(iD(q)));
    try
        Sp = conj_spectrum(z8, rv0(1:6), Tnd, cnd, mu, struct('K', K, 'nSub', nSub));
        R(q).nInterior   = Sp.nInterior;
        R(q).multiplicity = Sp.multiplicity;
        R(q).minRelSigma = Sp.minRel;
        R(q).nInteriorCand = Sp.nInteriorCand;
        R(q).nNearMiss   = Sp.nNearMiss;
        R(q).nZero       = Sp.nZero;
        R(q).nUnresolved = Sp.nUnresolved;     % located minimum inside the floor band
        R(q).nEndCand    = Sp.nEndCand;        % endpoint clusters, now refined too
        R(q).conjClear   = Sp.clear;
        R(q).candidates  = Sp.candidates;
        % the rank statement as a MEASURED margin: build C twice and let
        % Eckart-Young decide, instead of counting against a threshold.
        % H6 comes from the gates so its clearance is judged against this
        % arc's own Hamiltonian residual.
        gA = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, mu, ...
                                      struct('keepC', true, 'relTol', relPair(1)));
        gB = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, mu, ...
                                      struct('keepC', true, 'relTol', relPair(2)));
        Mg = lift_margin(gA.C, gB.C, z8(1:7), struct());
        R(q).h6margin    = gA.h6Margin;
        R(q).h6ok        = gA.h6Ok;
        R(q).h6clearance = gA.h6Clearance;
        R(q).liftMargin  = Mg.margin;
        R(q).liftCertified = Mg.certified;
        R(q).relTolPair  = relPair;
        R(q).done = true;
        intStr = '';
        for ci = find(~strcmp({Sp.candidates.class}, 'start'))
            intStr = [intStr sprintf(' [%s %s at t/tf %.4f, min %.1e x med]', ...
                      Sp.candidates(ci).class, Sp.candidates(ci).kind, ...
                      Sp.candidates(ci).tMinOverTf, Sp.candidates(ci).sigMinRel)]; %#ok<AGROW>
        end
        lg('  (%2d,%2d) interior %d, cand %d+%d (int+end): %d near-miss / %d zero / %d UNRESOLVED%s, minRel %.2e, H6 %.1fx %s, lift %.0fx %s', ...
           iD(q), iA(q), Sp.nInterior, Sp.nInteriorCand, Sp.nEndCand, Sp.nNearMiss, Sp.nZero, Sp.nUnresolved, intStr, Sp.minRel, ...
           gA.h6Margin, tern(gA.h6Ok, 'PASS', 'FAIL'), Mg.margin, tern(Mg.certified, 'CERT', 'uncert'));
    catch ME
        lg('  (%2d,%2d) THREW: %s', iD(q), iA(q), ME.message);
    end
    save(sideMat, 'R');                 % after EVERY entry, not at the end
    nThis = nThis + 1;
end

dn = [R.done];
S = struct('done', all(dn), 'nDone', nnz(dn), 'nTodo', nnz(~dn), 'rows', R, ...
           'worstSpectrum', maxOf([R(dn).nInterior]), ...
           'worstH6', minOf([R(dn).h6margin]), 'worstLift', minOf([R(dn).liftMargin]));
lg('second_order_pass: %d of %d measured; worst interior-crossing count %g, worst H6 margin %.2fx', ...
   S.nDone, n, S.worstSpectrum, S.worstH6);

if d('writeback', false)
    assert(S.done, 'census incomplete (%d todo) -- not writing back', S.nTodo);
    bak = [catMat '.bak_2nd'];
    if ~isfile(bak), copyfile(catMat, bak); end
    G = nan(size(s.has_solution));  M = G;  Hm = G;  Ci = G;  Cn = G;  Cz = G;  Cu = G;
    for q = 1:n
        G(iD(q), iA(q), iR(q))  = R(q).nInterior;
        M(iD(q), iA(q), iR(q))  = R(q).multiplicity;
        Hm(iD(q), iA(q), iR(q)) = R(q).h6margin;
        Ci(iD(q), iA(q), iR(q)) = R(q).nInteriorCand;
        Cn(iD(q), iA(q), iR(q)) = R(q).nNearMiss;
        Cz(iD(q), iA(q), iR(q)) = R(q).nZero;
        if isfield(R, 'nUnresolved') && ~isempty(R(q).nUnresolved)
            Cu(iD(q), iA(q), iR(q)) = R(q).nUnresolved;
        end
    end
    cat_.sheets(1).conj_interior = G;
    cat_.sheets(1).conj_multiplicity = M;
    cat_.sheets(1).conj_interior_cand = Ci;
    cat_.sheets(1).conj_near_miss = Cn;
    cat_.sheets(1).conj_zero = Cz;
    cat_.sheets(1).conj_unresolved = Cu;      % NaN = swept before 2026-09-11 (no such class)
    cat_.sheets(1).h6_margin = Hm;
    Lm = nan(size(s.has_solution));
    for q = 1:n, Lm(iD(q), iA(q), iR(q)) = R(q).liftMargin; end
    cat_.sheets(1).lift_margin = Lm;
    cat_.second_order = struct('date', datestr(now, 'yyyy-mm-dd'), ...
        'instruments', 'costate_common/{conj_spectrum,h6_margin,lift_margin}', ...
        'nSub', nSub, 'K', K, 'relTolPair', relPair, ...
        'meaning', ['conj_interior: interior sign changes of the quotiented determinant ' ...
                    'on a dense scan (K*nSub samples), 0 = none found; conj_multiplicity: ' ...
                    'INTERIOR candidates where two or more singular values collapsed together; ' ...
                    'conj_interior_cand / conj_near_miss / conj_zero / conj_unresolved: sigma_6 dips ' ...
                    'outside the start-up transient (endpoint dips INCLUDED since 2026-09-11), each LOCATED ' ...
                    'by shifted-grid refinement + golden section and judged against a numerical floor: ' ...
                    'zero (at the floor or a sign change), near-miss (positive minimum located), ' ...
                    'unresolved (inside the floor band; blocks a PASS). Entries swept before 2026-09-11 ' ...
                    'carry the OLD plateau classification and NaN in conj_unresolved: re-sweep to update; ' ...
                    'h6_margin: (c/T)/lambda_m(0), > 1 excludes ' ...
                    'the reduced problem''s spurious-zero mechanism; lift_margin: ' ...
                    'sigma_6 over the MEASURED error in the lift-space constraint matrix ' ...
                    '(Eckart-Young), > 1 certifies dim S = 1 rather than asserting it.']);
    Lout = struct(fn{1}, cat_);  save(catMat, '-struct', 'Lout');
    lg('[writeback] second-order measurements stored in %s (backup %s)', catMat, bak);
end
end

function legacyMatches(R, s, keys, sideMat)
% LEGACYMATCHES  Adopt an unkeyed sidecar only if it is provably the one
% the catalog was written from: every done record, and every written-back
% value equal. Errors second_order_pass:legacyMismatch otherwise.
% INPUTS: R sidecar records; s catalog sheet; keys [n x 3]; sideMat char.
pairs = {'h6_margin', 'h6margin'; 'lift_margin', 'liftMargin'; 'conj_interior', 'nInterior'; ...
         'conj_interior_cand', 'nInteriorCand'; 'conj_near_miss', 'nNearMiss'; 'conj_zero', 'nZero'};
if ~all([R.done])
    error('second_order_pass:legacyMismatch', 'legacy sidecar %s is incomplete: nothing to match against', sideMat);
end
for k = 1:size(pairs, 1)
    if ~isfield(s, pairs{k,1}) || ~isfield(R, pairs{k,2})
        error('second_order_pass:legacyMismatch', ...
              'legacy sidecar %s cannot be matched: %s is missing from the %s', sideMat, ...
              tern(isfield(s, pairs{k,1}), pairs{k,2}, pairs{k,1}), tern(isfield(s, pairs{k,1}), 'sidecar', 'catalog'));
    end
    for q = 1:size(keys, 1)
        a = s.(pairs{k,1})(keys(q,1), keys(q,2), keys(q,3));
        b = R(q).(pairs{k,2});
        if ~(isscalar(b) && isequaln(double(a), double(b)))
            error('second_order_pass:legacyMismatch', ...
                  'legacy sidecar %s disagrees with the catalog at record %d, %s', sideMat, q, pairs{k,1});
        end
    end
end
end

function v = maxOf(x), if isempty(x), v = NaN; else, v = max(x); end, end
function v = minOf(x), x = x(isfinite(x)); if isempty(x), v = NaN; else, v = min(x); end, end
function s = tern(c, a, b), if c, s = a; else, s = b; end, end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function logmsg(f, s)
% LOGMSG  Print, and append to a log file when one is named.  INPUTS: f; s.
fprintf('%s\n', s);
if ~isempty(f), fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
