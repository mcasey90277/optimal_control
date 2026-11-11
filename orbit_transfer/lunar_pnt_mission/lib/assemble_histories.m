function H = assemble_histories(M, P, assign, opts)
%% Purpose:
%
%   Puts the six spacecraft on ONE clock and ONE uniform time grid. Each
%   spacecraft coasts in its DRO until its departure, flies its polished
%   low-thrust transfer, then sits on its tulip slot. The hand-over
%   position continuity is MEASURED, not assumed (H.jumpKm).
%
%   Sampling rule (clock time t, ND):
%   (the third dimension of H.r and the columns of H.phase index the
%   SPACECRAFT 1..6, not the slot: spacecraft k sits on slot H.assign(k).)
%     t <  tDep   DRO coast   M.depState(mod(droPhase0(k) + t/tauD, 1))
%     tDep <= t < tArr   transfer   interp1(tDep + flight.t, flight.Y(:,1:3),
%                                   t, 'pchip'). Repeated flight times are
%                                   de-duplicated first (first kept); the
%                                   result must be strictly increasing
%                                   (assemble_histories:flightGrid). The
%                                   query is clamped to the last flight
%                                   sample, which is only ever a gap below
%                                   opts.epochTol: a larger mismatch
%                                   |tDep + flight.t(end) - tArr| throws
%                                   assemble_histories:epochs.
%     t >= tArr   on station  M.arrState(mod(slotPhase0(assign(k)) + t/tauT, 1))
%   M.depState / M.arrState (phase_state closures) accept a vector of
%   phases and return [6 x m]; the history uses them vectorised.
%
%   A spacecraft whose P.ok is false and which is not marked onStation is
%   REFUSED (assemble_histories:notAccepted): an unaccepted transfer is
%   never flown silently.
%
%   Modes: onSlots = true puts all six on their slots from t = 0 (no P; the
%   validation-gate mode; tEnd required). opts.onStation(k) = true puts
%   only spacecraft k on its slot from t = 0 (its P entry may be empty).
%
%% Inputs:
%
%  M                        struct                  mission_inputs() output
%  P                        [1 x 6] cell OR [1 x 6] struct array, or []
%                                                   polish_rendezvous
%                                                   results (.ok .tDep .tArr
%                                                   .flight.t [n x 1] ND from
%                                                   0 at departure,
%                                                   .flight.Y [n x 14]).
%                                                   Empty entries allowed
%                                                   where onStation(k).
%  assign                   [1 x 6]                 Slot of each spacecraft
%  opts                     struct (optional)
%    .dtSec                 scalar [300]            Grid step [s]
%    .serviceDays           scalar [30]             Days after the last
%                                                   arrival (default tEnd)
%    .tEnd                  scalar [ND]             Last epoch (default
%                                                   tLast + service); must
%                                                   be a positive scalar and
%                                                   >= tLast
%                                                   (assemble_histories:tEnd)
%    .onStation             [1 x 6] logical [false] On slot from t = 0
%    .onSlots               logical [false]         All six on slots
%    .epochTol              scalar [1e-8]           Largest tolerated
%                                                   |tDep + flight.t(end) -
%                                                   tArr| [ND]; larger throws
%                                                   (inf disables the check)
%
%% Outputs:
%
%  H.t          [N x 1]        ND, uniform from 0, step dtSec/tStar
%  H.r          [N x 3 x 6]    ND barycentric rotating positions; page k =
%                              spacecraft k
%  H.phase      [N x 6] int8   column k = spacecraft k; 1 DRO coast, 2 transfer, 3 on station
%  H.tDep       [1 x 6]        Departure epoch [ND]; NaN for spacecraft on
%                              station from t = 0 (0 is a legal departure)
%  H.tArr       [1 x 6]        Arrival epoch [ND]; NaN likewise
%  H.tFirst     scalar         Assembly window start: earliest tArr over the
%                              spacecraft that fly a transfer [ND] (NaN
%                              ignored); 0 if none
%  H.tLast      scalar         Assembly window end: latest tArr over the
%                              spacecraft that fly a transfer [ND]; 0 if none
%                              (the service window is [tLast, tLast +
%                              serviceDays]; grid ends are H.t(1), H.t(end))
%  H.tEnd       scalar         Requested tEnd [ND]
%  H.assign     [1 x 6]        Slot of each spacecraft
%  H.jumpKm     struct         .dep [1 x 6]: DRO state at tDep vs the
%                              flight's first sample; .arr [1 x 6]: the
%                              flight's last sample vs the slot state at
%                              tArr; position difference [km]. NaN for
%                              spacecraft on station from t = 0.
%  H.dtSec      scalar         Grid step [s]
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
    M0 = mission_inputs();
    H = assemble_histories(M0, [], 1:6, struct('onSlots', true, ...
                           'tEnd', 5*86400/M0.tStar));
    fprintf('on-slots demo: N = %d, r(1,:,1) = [%.4f %.4f %.4f] ND\n', ...
            numel(H.t), H.r(1, :, 1));
    return;
end

if nargin < 4, opts = struct(); end
dtSec        = getopt(opts, 'dtSec', 300);
serviceDays  = getopt(opts, 'serviceDays', 30);
onSlots      = getopt(opts, 'onSlots', false);
epochTol     = getopt(opts, 'epochTol', 1e-8);
onStation    = logical(getopt(opts, 'onStation', false(1, 6)));
if onSlots, onStation = true(1, 6); end
assign = assign(:).';
if numel(assign) ~= 6 || numel(onStation) ~= 6
    error('assemble_histories:shape', 'assign and onStation must have 6 entries');
end

nSc = 6;
tDep = nan(1, nSc);  tArr = nan(1, nSc);
Pk   = cell(1, nSc);
tFl  = cell(1, nSc);  YFl = cell(1, nSc);
for kk = 1:nSc
    if onStation(kk), continue; end
    if isstruct(P) && numel(P) == nSc
        Pk{kk} = P(kk);
    elseif iscell(P) && numel(P) == nSc
        Pk{kk} = P{kk};
    else
        error('assemble_histories:shape', ...
              'P must be a [1 x 6] cell or struct array (spacecraft %d)', kk);
    end
    if isempty(Pk{kk}) || ~isfield(Pk{kk}, 'ok') || ~Pk{kk}.ok
        error('assemble_histories:notAccepted', ...
              'spacecraft %d has no accepted transfer and is not onStation', kk);
    end
    tDep(kk) = Pk{kk}.tDep;
    tArr(kk) = Pk{kk}.tArr;
    if tDep(kk) < 0 || tArr(kk) <= tDep(kk)
        error('assemble_histories:epochs', ...
              'spacecraft %d: tDep %.6g, tArr %.6g', kk, tDep(kk), tArr(kk));
    end
    tf0 = Pk{kk}.flight.t(:);
    Yf0 = Pk{kk}.flight.Y(:, 1:3);
    keep = [true; diff(tf0) ~= 0];
    tf0 = tf0(keep);  Yf0 = Yf0(keep, :);
    if numel(tf0) < 2 || any(diff(tf0) <= 0)
        error('assemble_histories:flightGrid', ...
              'spacecraft %d: flight times are not strictly increasing', kk);
    end
    tFl{kk} = Pk{kk}.tDep + tf0;
    YFl{kk} = Yf0;
    gap = abs(tFl{kk}(end) - tArr(kk));
    if gap > epochTol
        error('assemble_histories:epochs', ...
              ['spacecraft %d: flight ends %.3g ND from tArr (tolerance %.3g): ' ...
               'bad tArr or flight length'], kk, gap, epochTol);
    end
end

if any(~onStation)
    tFirst = min(tArr(~onStation));  tLast = max(tArr(~onStation));
else
    tFirst = 0;  tLast = 0;
end
dt = dtSec/M.tStar;
if isfield(opts, 'tEnd') && ~isempty(opts.tEnd)
    tEnd = opts.tEnd;
elseif onSlots
    error('assemble_histories:tEnd', 'onSlots mode needs opts.tEnd');
else
    tEnd = tLast + serviceDays*86400/M.tStar;
end
if ~isscalar(tEnd) || ~isfinite(tEnd) || tEnd <= 0 || tEnd < tLast
    error('assemble_histories:tEnd', ...
          'tEnd must be a positive scalar >= tLast = %.6g ND', tLast);
end
nn = floor(tEnd/dt + 1e-9) + 1;
t  = (0:nn-1).'*dt;

r     = zeros(nn, 3, nSc);
phase = zeros(nn, nSc, 'int8');
jDep  = nan(1, nSc);
jArr  = nan(1, nSc);
for kk = 1:nSc
    sSlot = @(tt) mod(M.slotPhase0(assign(kk)) + tt/M.tauT, 1);
    if onStation(kk)
        rv = M.arrState(sSlot(t));
        r(:, :, kk) = rv(1:3, :).';
        phase(:, kk) = 3;
        continue;
    end
    Pc  = Pk{kk};
    tF  = tFl{kk};
    Yf  = YFl{kk};
    iC  = t < Pc.tDep;
    iT  = t >= Pc.tDep & t < Pc.tArr;
    iS  = t >= Pc.tArr;
    sDro = @(tt) mod(M.droPhase0(kk) + tt/M.tauD, 1);
    if any(iC)
        rv = M.depState(sDro(t(iC)));
        r(iC, :, kk) = rv(1:3, :).';
    end
    if any(iT)
        tq = min(t(iT), tF(end));
        r(iT, :, kk) = interp1(tF, Yf, tq, 'pchip');
    end
    if any(iS)
        rv = M.arrState(sSlot(t(iS)));
        r(iS, :, kk) = rv(1:3, :).';
    end
    phase(iC, kk) = 1;  phase(iT, kk) = 2;  phase(iS, kk) = 3;
    % hand-over continuity [km]
    rd = M.depState(sDro(Pc.tDep));
    jDep(kk) = sqrt(sum((rd(1:3) - Yf(1, :).').^2))*M.lStar;
    ra = M.arrState(sSlot(Pc.tArr));
    jArr(kk) = sqrt(sum((ra(1:3) - Yf(end, :).').^2))*M.lStar;
end

H.t      = t;
H.r      = r;
H.phase  = phase;
H.tDep   = tDep;
H.tArr   = tArr;
H.tFirst = tFirst;
H.tLast  = tLast;
H.tEnd   = tEnd;
H.assign = assign;
H.jumpKm = struct('dep', jDep, 'arr', jArr);
H.dtSec  = dtSec;
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
