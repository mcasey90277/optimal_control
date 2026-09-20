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
%   .sD [] an explicit list of departure phases (any spacing); else the
%   lattice .sDorigin [0] + (0:nD-1)/nD, .nD [12]
%   .sD0 [0] .thrustN [0.070] .ispS [900] .m0kg [150]
%   .tauDRO [1] .NpTulip [7] .pmTulip [-1]
%   .families [] a family_map output: every entry is then stamped with
%   the code of the extremal family it belongs to (Q.FAM), and the map
%   itself ships with the sheet (Q.families)
%
%% Outputs:
%
%  Q                        struct                  .OK [nD x nA x 1]
%                                                   .Z8 [8 x nD x nA x 1]
%                                                   .TF [nD x nA x 1] ND
%                                                   .CONJ (int8, -1 = no
%                                                   entry) .MINLV .MINQ
%                                                   .DIMS the verdicts
%                                                   .FAM (int8 family codes,
%                                                   with .families) when a
%                                                   map was given
%                                                   .sD .sA .rungs .meta
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
nA = numel(S.sA);

% THE PROBLEM IDENTITY IS AUTHORITATIVE. Engine, orbits and the departure
% phase come from what was CERTIFIED, never from fresh option defaults --
% otherwise packaging can label a real trajectory with the wrong physics,
% and a sheet certified at one departure phase lands on another's row.
% Options may ASSERT a value; disagreement is an error, not an override.
% (Astra chain review 2026-09-10.)
assert(isfield(S, 'problem') && isstruct(S.problem), ...
       ['sheet carries no problem identity: rebuild it with build_arrival_sheet, ' ...
        'which stamps S.problem from the certified setup']);
P = S.problem;
tStar = P.tStar;  lStar = P.lStar;
nD = d('nD', 12);
assertMatch(opts, P, {'thrustN', 'ispS', 'm0kg', 'tauDRO', 'NpTulip', 'pmTulip'});
if isfield(opts, 'sD0'), assertNear(opts.sD0, P.sD, 'sD0'); end
sD0 = P.sD;                                   % the CERTIFIED departure phase

Q = struct();
% THE GRID ORIGIN IS NOT THE CERTIFIED PHASE. Building the departure grid
% from sD0 made idxOf(Q.sD, sD0) select row 1 by construction, so a sheet
% certified at any departure phase was placed at departure zero. The grid is
% the canonical k/nD one; the certified phase must LAND on it.
if isfield(opts, 'sD') && ~isempty(opts.sD)
    Q.sD = mod(opts.sD(:).', 1);  nD = numel(Q.sD);       % the caller's list
else
    Q.sD = mod(d('sDorigin', 0) + (0:nD-1)/nD, 1);
end
Q.sA = S.sA(:)';
Q.rungs = P.thrustN;
Q.OK = false(nD, nA, 1);
Q.TF = nan(nD, nA, 1);
Q.Z8 = nan(8, nD, nA, 1);
% THE VERDICTS TRAVEL WITH THE ENTRY. Certification already establishes the
% conjugate test and the three sufficiency-hypothesis gates for every point
% here; discarding them at packaging made the catalog report a conjugate
% census of 0/0/0 on entries that were all conjugate-certified, so a
% recipient reading conj_pass would have seen nothing. -1 = no entry (the
% same convention conj_catalog_pass uses).
Q.CONJ  = -ones(nD, nA, 1, 'int8');
Q.MINLV = nan(nD, nA, 1);
Q.MINQ  = nan(nD, nA, 1);
Q.DIMS  = nan(nD, nA, 1);
% the junction count the conjugate test SAMPLED AT is part of what the
% verdict means (it is a sign test at K-1 interior junctions), so it travels
% with the verdict rather than being reconstructed later
Q.KJ    = nan(nD, nA, 1);
% THE FAMILY of every entry (family_map codes: 1..nFam mapped, -1 a root no
% arc passes through, -2 a rib whose spine root is unidentified, 0 none)
F = d('families', []);
if ~isempty(F)
    Q.FAM = zeros(nD, nA, 1, 'int8');
    Q.families = rmfield(F, intersect(fieldnames(F), {'attach', 'ribFamily'}));
    % WHAT family_index IS, said in the file that carries it
    Q.families.note = ['family_index is PROVENANCE: which continuation arc (family) an entry was found on, ' ...
        'as family_map attached it. Spine entries are attached by flight time and costates; rib entries by ' ...
        'flight time alone (code -2 = unidentified), which is NOT reproducible between builds (FINDINGS 82). ' ...
        'Nothing in the pipeline reads it. To decide whether two neighbouring entries may be interpolated ' ...
        'between, use the branch map of phase_transversality_check (.safeD .safeA .branch), which is computed ' ...
        'from the entries themselves.'];
end
% THE ENTRY NOTES: how each entry was found and what the certifier remarked,
% one short text per cell (the family's label in front when a map exists)
Q.NOTE = repmat({''}, nD, nA);

% ---- the spine: the certified minimum at each arrival phase, at sD0 -----
iD0 = idxOf(Q.sD, sD0);
assert(~isempty(iD0), 'the certified departure phase %.6f is not on the %d-point grid', sD0, nD);
for j = 1:nA
    if ~isfinite(S.TF(j)), continue, end
    % EXPORT FROM THE CERTIFICATE, not from the summary. A finite S.TF is a
    % summary claim; without the certificate behind it there is nothing to
    % ship, and t_f is taken from the certificate's own z(8) so the two
    % representations cannot disagree.
    c = S.cand{j};
    if isempty(c), continue, end
    k = find([c.ok] & abs([c.tfDays] - S.TF(j)) < 1e-9, 1);
    if isempty(k), continue, end
    if ~usableEntry(c(k)), continue, end
    Q.OK(iD0, j, 1) = true;
    Q.TF(iD0, j, 1) = c(k).z(8);
    Q.Z8(:, iD0, j, 1) = c(k).z(:);
    Q = putVerdicts(Q, iD0, j, c(k));
    code = 0;
    if ~isempty(F), Q.FAM(iD0, j, 1) = int8(F.columns(j).family);  code = Q.FAM(iD0, j, 1); end
    Q.NOTE{iD0, j} = noteOf(c(k), F, code);
end

% ---- the ribs: certified departure points off the spine ----------------
if nargin >= 2 && ~isempty(ribs)
    if ~iscell(ribs), ribs = {ribs}; end
    for k = 1:numel(ribs)
        R = ribs{k};
        if ~isfield(R, 'pts') || isempty(R.pts), continue, end
        ribCode = int8(0);
        if ~isempty(F)
            okp = R.pts(arrayfun(@usableEntry, R.pts));
            tfp = [okp(1:min(3, end)).tfDays];   % the points nearest the spine
            ribCode = int8(F.ribFamily(R.sA, tfp));
        end
        for m = 1:numel(R.pts)
            Pt = R.pts(m);
            if ~usableEntry(Pt), continue, end
            iD = idxOf(Q.sD, Pt.sD);  iA = idxOf(Q.sA, Pt.sA);
            assert(~isempty(iD) && ~isempty(iA), 'rib point (%.4f, %.4f) is off the grid', Pt.sD, Pt.sA);
            % keep the faster. Compared only AFTER usableEntry, because
            % `existingTF <= NaN` is false and an unusable point would
            % otherwise displace a good one.
            if Q.OK(iD, iA, 1) && Q.TF(iD, iA, 1) <= Pt.z(8), continue, end
            Q.OK(iD, iA, 1) = true;
            Q.TF(iD, iA, 1) = Pt.z(8);
            Q.Z8(:, iD, iA, 1) = Pt.z(:);
            Q = putVerdicts(Q, iD, iA, Pt);
            if ~isempty(F), Q.FAM(iD, iA, 1) = ribCode; end
            Q.NOTE{iD, iA} = noteOf(Pt, F, ribCode);
        end
    end
end

% ---- meta: what the packager reads -------------------------------------
% THE PERIOD COMES FROM THE SHEET, not from a literal. Both fields used to
% be the hardcoded 5*2*pi/6 beside a P.NpTulip read straight from the
% certified problem, so a catalog built at another petal count would have
% carried a 7-petal period as its own label. Identical for Np = 7.
tauT = P.periodTulip;
Q.meta = struct('muStar', P.muStar, 'lStar', lStar, 'tStar', tStar, ...
    'ispS', P.ispS, 'm0kg', P.m0kg, 'tauDRO', P.tauDRO, ...
    'depFamily', 'dro', 'depParams', struct('tau', P.tauDRO), ...
    'arrFamily', 'tulip', 'arrParams', struct('Np', P.NpTulip, 'pm', P.pmTulip, 'tau', tauT), ...
    'NpTulip', P.NpTulip, 'pmTulip', P.pmTulip, 'periodTulip', tauT);
Q.problem = P;                                 % identity ships with the sheet

if ~isempty(outMat), save(outMat, '-struct', 'Q'); end
end

function ok = usableEntry(C)
% USABLEENTRY  A certificate is exportable only if it says ok AND carries
% numbers worth shipping: a real finite 8-vector with a positive final time.
% `P.ok = true` with `z(8) = NaN` used to pass, and could displace a good
% entry because `existingTF <= NaN` is false.  INPUTS: C.  OUTPUTS: ok.
ok = false;
if ~isstruct(C) || ~isscalar(C) || ~isfield(C, 'ok') || ~isfield(C, 'z'), return, end
[okf, f] = scalar_verdict(C.ok);
if ~okf || f ~= 1, return, end
z = C.z;
if ~isnumeric(z) || numel(z) ~= 8 || ~all(isfinite(z(:))) || ~isreal(z), return, end
if ~(z(8) > 0), return, end
ok = true;
end

function assertMatch(opts, P, names)
% ASSERTMATCH  Options may assert an identity value, never replace it.
% INPUTS: opts; P; names.
for k = 1:numel(names)
    if isfield(opts, names{k}) && ~isempty(opts.(names{k}))
        assertNear(opts.(names{k}), P.(names{k}), names{k});
    end
end
end

function assertNear(a, b, name)
% ASSERTNEAR  INPUTS: a; b; name.
assert(abs(a - b) <= 1e-12*max(abs(b), 1), ...
    ['opts.%s = %g contradicts the CERTIFIED problem identity (%g). ' ...
     'Packaging may assert the identity, not change it.'], name, a, b);
end

function s = noteOf(C, F, code)
% NOTEOF  An entry's note: the family label (when a map exists), then the
% producer's and certifier's text.  INPUTS: C; F ([] = no map); code.
% OUTPUTS: s char.
parts = {};
if ~isempty(F)
    if code >= 1, parts{end+1} = ['family ' F.families(code).label];
    elseif code == -1, parts{end+1} = 'family: unattached root (no arc passes through it)';
    elseif code == -2, parts{end+1} = 'family: spine unidentified';
    end
end
if isfield(C, 'note') && ~isempty(C.note), parts{end+1} = C.note; end
s = strjoin(parts, ' | ');
end

function Q = putVerdicts(Q, iD, iA, C)
% PUTVERDICTS  Store one entry's conjugate verdict and hypothesis gates.
% INPUTS: Q; iD; iA; C (a certify_root output).  OUTPUTS: Q.
if isfield(C, 'conj') && ~isempty(C.conj), Q.CONJ(iD, iA, 1) = int8(C.conj); end
if isfield(C, 'Y') && ~isempty(C.Y), Q.KJ(iD, iA, 1) = size(C.Y, 2); end
if isfield(C, 'g') && isstruct(C.g) && ~isempty(C.g)
    Q.MINLV(iD, iA, 1) = C.g.minLamV;
    Q.MINQ(iD, iA, 1)  = C.g.minQmt;
    Q.DIMS(iD, iA, 1)  = C.g.dimS;
end
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
