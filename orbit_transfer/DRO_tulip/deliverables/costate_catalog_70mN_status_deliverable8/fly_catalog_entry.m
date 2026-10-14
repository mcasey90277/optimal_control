function F = fly_catalog_entry(cat, entry, opts)
%% Purpose:
%
%   Fly one library entry and hand back everything the summary printer and
%   the 3D figure need: from its multiple-shooting junction states when it
%   has them (accurate over 17-26 days), otherwise from its z8 in one
%   flight. The arrival target is the catalog's own endpoint rule
%   (catalog_endpoints) at the entry's phase fractions.
%
%  ASSUMPTIONS / NOTES:
%
% • Phases are FRACTIONS of a period and wrap (catalog_endpoints).
% • Delta-V is the rocket equation cnd*log(1/mf), as fly_transfer forms it.
% • Thrust direction u = -lambda_v/|lambda_v| (all-burn min-time).
% • entry.status / .statusName / .reason are carried from a transfers_at
%   row; absent, status is NaN, statusName and reason 'not annotated'.
%
%% Inputs:
%
%  cat                      struct                  library catalog
%  entry                    struct                  a transfers_at row, OR
%                                                   struct('iD',,'iA',) (the
%                                                   primary of that cell of
%                                                   sheet 1), OR struct('z8',
%                                                   ,'sD',,'sA',)
%  opts                     struct (optional)       .method ['auto'] |
%                                                   'junctions' | 'z8';
%                                                   .nSample [50] per
%                                                   segment; .B [] a
%                                                   catalog_endpoints result
%
%% Outputs:
%
%  F                        struct                  .method .reasonMethod
%                                                   .sD .sA .z .z8 .t .tDays
%                                                   .Y [n x 14] .u [n x 3]
%                                                   .massKg .tfNd .tfDays
%                                                   .dvKms .propellantKg
%                                                   .finalMassKg .flyKm
%                                                   .flyVms .defectKm
%                                                   .defectVms (empty for z8)
%                                                   .nSamples .admissibility
%                                                   .status .statusName
%                                                   .reason .rv0 .rvf .B
%
%% Revision History:
%  M. Casey                                                   (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3 || isempty(opts), opts = struct(); end
method  = 'auto';  if isfield(opts, 'method') && ~isempty(opts.method), method = opts.method; end
nSample = 50;      if isfield(opts, 'nSample') && ~isempty(opts.nSample), nSample = opts.nSample; end

%% Resolve the entry form:
if isfield(entry, 'iD') && isfield(entry, 'iA') && ~isfield(entry, 'z8')
    sh = cat.sheets(1);
    T  = transfers_at(cat, sh.sD_frac(entry.iD), sh.sA_frac(entry.iA));
    ip = find(strcmp({T.kind}, 'primary'), 1);
    if isempty(ip), error('fly_catalog_entry:noPrimary', 'cell (%d,%d) has no primary', entry.iD, entry.iA); end
    entry = T(ip);
end
if ~isfield(entry, 'z8') || numel(entry.z8) ~= 8
    error('fly_catalog_entry:entry', 'the entry needs a z8 of 8 elements');
end
z8 = entry.z8(:);
J  = [];  if isfield(entry, 'junctions'), J = entry.junctions; end
hasJ = ~isempty(J) && size(J, 1) == 14;
badShape = ~isempty(J) && size(J, 1) ~= 14;
if ~isfield(entry, 'sD') || ~isfield(entry, 'sA')
    error('fly_catalog_entry:entry', 'the entry needs phase fractions sD and sA');
end

switch method
    case 'auto'
        if hasJ, method = 'junctions';  reasonMethod = 'auto: junctions';
        elseif badShape
            method = 'z8';
            reasonMethod = sprintf('auto: junctions of unexpected size [%d x %d], flown from z8', size(J, 1), size(J, 2));
        else,    method = 'z8';         reasonMethod = 'auto: no junctions stored, flown from z8';
        end
    case 'junctions'
        if badShape
            error('fly_catalog_entry:junctionsShape', 'junctions must be 14 x K, got [%d x %d]', size(J, 1), size(J, 2));
        end
        if ~hasJ, error('fly_catalog_entry:noJunctions', 'the entry carries no junction states'); end
        reasonMethod = 'requested';
    case 'z8'
        reasonMethod = 'requested';
    otherwise
        error('fly_catalog_entry:method', 'unknown method %s', method);
end

%% Endpoints and physics:
if isfield(opts, 'B') && ~isempty(opts.B), B = opts.B; else, B = catalog_endpoints(cat); end
lStar = cat.constants.lStar_km;  tStar = cat.constants.tStar_s;  m0 = cat.thruster.m0_kg;
rv0 = B.stateD(entry.sD);  rv0 = rv0(1:6);
rvf = B.stateA(entry.sA);  rvf = rvf(1:6);
phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'mu', B.mu, 'muStar', B.mu, ...
              'lStar', lStar, 'tStar', tStar, 'm0kg', m0);

%% Fly:
if strcmp(method, 'junctions')
    P = entry_thrust_program(z8, J, rv0, rvf, phys, struct('nSample', nSample));
    t = P.t;  Y = P.Y;  u = P.u;  massKg = P.massKg;
    defectKm = P.defectKm;  defectVms = P.defectVms;
    flyKm = P.flyKm;  flyVms = P.flyVms;
else
    Ff = fly_transfer(z8, rv0, rvf, phys);
    t = Ff.t;  Y = Ff.Y;
    lv = Y(:, 11:13);
    u = -lv ./ sqrt(sum(lv.^2, 2));
    massKg = Y(:, 7)*m0;
    defectKm = [];  defectVms = [];
    flyKm = Ff.flyKm;  flyVms = Ff.flyVms;
end

%% Assemble:
mf = Y(end, 7);
F = struct('method', method, 'reasonMethod', reasonMethod, 'sD', entry.sD, 'sA', entry.sA, ...
           'z', z8, 'z8', z8, 't', t, 'tDays', t*tStar/86400, 'Y', Y, 'u', u, 'massKg', massKg, ...
           'tfNd', z8(8), 'tfDays', z8(8)*tStar/86400, ...
           'dvKms', B.cnd*log(1/mf)*lStar/tStar, 'propellantKg', (1 - mf)*m0, ...
           'finalMassKg', mf*m0, 'flyKm', flyKm, 'flyVms', flyVms, ...
           'defectKm', defectKm, 'defectVms', defectVms, 'nSamples', numel(t));
F.admissibility = validate_flight(t, Y, z8(8), B.Tnd, B.cnd, B.mu, lStar);
F.status = NaN;  F.statusName = 'not annotated';  F.reason = 'not annotated';
if isfield(entry, 'status') && ~isempty(entry.status) && isfinite(entry.status)
    F.status = entry.status;
    K = status_key();
    kk = find(K.codes == entry.status, 1);
    if ~isempty(kk), F.statusName = K.names{kk}; end
end
if isfield(entry, 'statusName') && ~isempty(entry.statusName), F.statusName = entry.statusName; end
if isfield(entry, 'reason'), F.reason = entry.reason; end
F.rv0 = rv0;  F.rvf = rvf;  F.B = B;
end
