function F = regime_features(A, rec, fam)
%% Purpose:
%
%   Feature vector of ONE continuation arm, for the REGIME MAP: what the
%   walk achieved, what it cost, and what the trajectory looked like where
%   the walk stopped -- the quantities the regime hypotheses are stated in:
%
%     H1  huber walls at a GRAZING-RISK switch-structure change
%         -> nCross, minAbsDQdt, nGraze, grazeGap
%     H2  eps fails at the floor on MANY-SWITCH cells (ramp slope 1/2p)
%         -> nCross, coastFrac, revs, tfDays, nFail
%     H3  huberc's early walls are unexplained
%         -> condJFloor, condJRatio, plateauFrac (fraction of the arc with
%            Q within the ramp of 1: the band the delta-ramp must spread)
%
%   Every family is put on ONE sharpness axis, the effective ramp width in
%   Q: eps 2p (its law ramps over [1-p, 1+p]), huberc delta, huber 0 (the
%   jump). "Reached the floor" then means the same thing across families.
%
%   Trajectory features are measured at the LAST ACCEPTED rung -- the state
%   from which the walk failed to continue, which is what a wall diagnosis
%   needs (a walled arm has no floor solution to look at).
%
%% Inputs:
%
%  A                        struct                  One arm of run_minfuel_race
%                                                   (.p .delta .mf .coastFrac
%                                                   .condJ .Y .nFail .nBisect
%                                                   .retired .wallTotal)
%
%  rec                      struct                  Its seed record (pilot
%                                                   format: .iD .iA .gam .tf
%                                                   .rv0 .rvf .Tmax .c .muStar)
%
%  fam                      char                    'eps' | 'huber' | 'huberc'
%
%% Outputs:
%
%  F                        struct                  .family .cellIdx .gamma
%                                                   .tfDays .outcome
%                                                   ('floor'|'wall'|'none')
%                                                   .pFloor .deltaFloor
%                                                   .rampWidth .mf .coastFrac
%                                                   .nRung .nFail .nBisect
%                                                   .wallMin .condJFloor
%                                                   .condJRatio .nCross
%                                                   .minAbsDQdt .nGraze
%                                                   .grazeGap .plateauFrac
%                                                   (|Q-1| <= 0.05 fraction)
%                                                   .revs .nCrossStart
%                                                   .dSwitch (switches the
%                                                   walk had to create)
%
%% Revision History:
%  M. Casey                                                   (c) 09/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

pFloorTarget = 1e-3;                                  % the schedule's last rung
floorTol     = 1.5;                                   % a bisected last rung lands
                                                      % just ABOVE the target (eps
                                                      % on (1,2)@1.223 stops at
                                                      % 0.001311 and is the shipped
                                                      % catalog entry); within this
                                                      % factor the walk has arrived

F = struct('family', fam, 'cellIdx', [rec.iD rec.iA], 'gamma', rec.gam, ...
    'tfDays', rec.tf*4.3481, 'outcome', 'none', 'pFloor', NaN, 'deltaFloor', NaN, ...
    'rampWidth', NaN, 'mf', NaN, 'coastFrac', NaN, 'nRung', 0, ...
    'nFail', A.nFail, 'nBisect', A.nBisect, 'wallMin', NaN, ...
    'condJFloor', NaN, 'condJRatio', NaN, 'nCross', NaN, 'minAbsDQdt', NaN, ...
    'nGraze', NaN, 'grazeGap', NaN, 'plateauFrac', NaN, 'revs', NaN, ...
    'nCrossStart', NaN, 'dSwitch', NaN);
if isfield(A, 'wallTotal'), F.wallMin = A.wallTotal/60; end
if isempty(A.p), return, end

F.nRung  = numel(A.p);
F.pFloor = A.p(end);
F.mf     = A.mf(end);
F.coastFrac = A.coastFrac(end);
if ~isempty(A.delta), F.deltaFloor = A.delta(end); end
if ~isempty(A.condJ)
    F.condJFloor = A.condJ(end);
    F.condJRatio = A.condJ(end) / A.condJ(1);
end

% --- the common sharpness axis (effective ramp width in Q) ---------------
switch fam
    case 'eps',    F.rampWidth = 2*F.pFloor;    F.deltaFloor = NaN;
    case 'huber',  F.rampWidth = 0;             F.deltaFloor = NaN;
    case 'huberc', F.rampWidth = F.deltaFloor;
    otherwise, error('regime_features:family', 'unknown family %s', fam);
end

% --- outcome: did the p-walk reach the schedule's last rung? -------------
% Outcome is read on p, the ONE continuation parameter all three families
% walk (ramp width cannot serve: huber's is 0 at every p by construction).
if F.pFloor <= pFloorTarget*floorTol, F.outcome = 'floor'; else, F.outcome = 'wall'; end

% --- trajectory features at the LAST ACCEPTED rung -----------------------
sm = struct('family', fam, 'p', F.pFloor);
if strcmp(fam, 'huberc'), sm.delta = F.deltaFloor; end
Y1 = A.Y{end}(:,1);
try
    D = huber_switch_diag(Y1, rec.tf, rec.Tmax, rec.c, rec.muStar, sm);
    F.nCross = D.nCross;
    F.minAbsDQdt = D.minAbsDQdt;
    F.nGraze = numel(D.graze);
    if F.nGraze > 0, F.grazeGap = min(abs([D.graze.Q] - 1)); end
    % plateau: fraction of the arc with Q within a FIXED band of 1 -- a
    % geometry feature of the cell, so it must not depend on the family's
    % own ramp width (which is the axis, not the measurement)
    F.plateauFrac = mean(abs(D.Q - 1) <= 0.05);
catch
    % a stalled arm can carry a state the propagator refuses; features stay NaN
end
F.revs = moonRevs(Y1, rec, sm);

% --- how many switches did the walk have to CREATE? ----------------------
% The energy seed (first accepted rung, p ~ 1) has its own switch structure;
% the bang-bang limit has another. The difference is the structural work the
% continuation must do, and it is the quantity the wall hypotheses turn on.
sm0 = struct('family', fam, 'p', A.p(1));
if strcmp(fam, 'huberc') && ~isnan(A.delta(1)), sm0.delta = A.delta(1); end
try
    D0 = huber_switch_diag(A.Y{1}(:,1), rec.tf, rec.Tmax, rec.c, rec.muStar, sm0);
    F.nCrossStart = D0.nCross;
    F.dSwitch = F.nCross - F.nCrossStart;
catch
end
end

% ------------------------------------------------------------------------
function n = moonRevs(Y1, rec, sm)
% MOONREVS  Revolutions of the arc about the Moon in the rotating frame:
% the total unwrapped swept angle of (r - r_Moon) divided by 2*pi, flown in
% the ARM's OWN smoothing family.
% INPUTS: Y1 [14x1] initial PMP state; rec (seed record); sm (smoothing).
% OUTPUTS: n.
n = NaN;
try
    [~, ~, T_, Yt] = cr3bp_minfuel_prop(rec.tf, Y1, false, rec.Tmax, rec.c, rec.muStar, sm);
    [~, iu] = unique(T_);  Yt = Yt(iu, :);
    dx = Yt(:,1) - (1 - rec.muStar);  dy = Yt(:,2);
    n = abs(unwrap(atan2(dy, dx)));
    n = (max(n) - min(n)) / (2*pi);
catch
end
end
