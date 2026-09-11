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
%
%% Outputs:
%
%  S                        struct                  .done .nDone .nTodo
%                                                   .worstSpectrum .worstH6
%                                                   .worstLift .rows
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
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
g0 = 9.80665*tStar^2/(1000*lStar);
cnd = (isp/tStar)*g0;
Tnd = (cat_.rungs_N(1)/m0)*tStar^2/(lStar*1000);

ob = struct('muStar', mu, 'lStar', lStar, 'tStar', tStar, 'tauDRO', s.tauDRO, ...
            'NpTulip', s.Np, 'tauTulip', s.period_tulip_nd, 'pmTulip', s.pm, ...
            'ispS', isp, 'm0kg', m0);
[tD, rvD, ~, ~] = ladder_endpoints(ob);
stD = @(x) interp1(tD, rvD, mod(x,1)*tD(end), 'spline')';

[iD, iA, iR] = ind2sub(size(s.has_solution), find(s.has_solution));
n = numel(iD);
% resume from the sidecar
if isfile(sideMat)
    P = load(sideMat);  R = P.R;
    lg('resuming from %s (%d of %d already measured)', sideMat, nnz([R.done]), n);
else
    R = struct('done', num2cell(false(1,n)), 'nInterior', [], 'multiplicity', [], ...
               'minRelSigma', [], 'h6margin', [], 'h6ok', [], 'liftMargin', [], ...
               'liftCertified', []);
    R = R(:).';
end

nThis = 0;
for q = 1:n
    if R(q).done, continue, end
    if toc(t0) > batchSec, lg('batch budget reached, exiting cleanly'); break, end
    if nThis >= maxEnt, lg('maxEntries reached'); break, end
    z8 = s.z8(:, s.entry_index(iD(q), iA(q), iR(q)));
    rv0 = stD(s.sD_frac(iD(q)));
    try
        Sp = conj_spectrum(z8, rv0(1:6), Tnd, cnd, mu, struct('K', K, 'nSub', nSub));
        H6 = h6_margin(z8, Tnd, cnd);
        R(q).nInterior   = Sp.nInterior;
        R(q).multiplicity = Sp.multiplicity;
        R(q).minRelSigma = Sp.minRel;
        R(q).h6margin    = H6.margin;
        R(q).h6ok        = H6.ok;
        % the rank statement as a MEASURED margin: build C twice and let
        % Eckart-Young decide, instead of counting against a threshold
        gA = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, mu, ...
                                      struct('keepC', true));
        gB = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, mu, ...
                                      struct('keepC', true, 'relTol', 1e-7));
        Mg = lift_margin(gA.C, gB.C, z8(1:7), struct());
        R(q).liftMargin  = Mg.margin;
        R(q).liftCertified = Mg.certified;
        R(q).done = true;
        lg('  (%2d,%2d) interior %d, mult %d, minRel %.2e, H6 %.1fx %s, lift %.0fx %s', ...
           iD(q), iA(q), Sp.nInterior, Sp.multiplicity, Sp.minRel, H6.margin, ...
           tern(H6.ok, 'PASS', 'FAIL'), Mg.margin, tern(Mg.certified, 'CERT', 'uncert'));
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
    G = nan(size(s.has_solution));  M = G;  Hm = G;
    for q = 1:n
        G(iD(q), iA(q), iR(q))  = R(q).nInterior;
        M(iD(q), iA(q), iR(q))  = R(q).multiplicity;
        Hm(iD(q), iA(q), iR(q)) = R(q).h6margin;
    end
    cat_.sheets(1).conj_interior = G;
    cat_.sheets(1).conj_multiplicity = M;
    cat_.sheets(1).h6_margin = Hm;
    Lm = nan(size(s.has_solution));
    for q = 1:n, Lm(iD(q), iA(q), iR(q)) = R(q).liftMargin; end
    cat_.sheets(1).lift_margin = Lm;
    cat_.second_order = struct('date', datestr(now, 'yyyy-mm-dd'), ...
        'instruments', 'costate_common/{conj_spectrum,h6_margin,lift_margin}', ...
        'nSub', nSub, 'K', K, ...
        'meaning', ['conj_interior: interior sign changes of the quotiented determinant ' ...
                    'on a dense scan (K*nSub samples), 0 = none found; conj_multiplicity: ' ...
                    'candidates where two or more singular values collapsed together, which ' ...
                    'a determinant cannot see; h6_margin: (c/T)/lambda_m(0), > 1 excludes ' ...
                    'the reduced problem''s spurious-zero mechanism; lift_margin: ' ...
                    'sigma_6 over the MEASURED error in the lift-space constraint matrix ' ...
                    '(Eckart-Young), > 1 certifies dim S = 1 rather than asserting it.']);
    Lout = struct(fn{1}, cat_);  save(catMat, '-struct', 'Lout');
    lg('[writeback] second-order measurements stored in %s (backup %s)', catMat, bak);
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
