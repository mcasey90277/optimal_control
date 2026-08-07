function [lam, tStations, diag_] = duals_to_costates(spec)
%% Purpose:
%
%   Maps the KKT multipliers of a direct transcription's DEFECT constraints
%   to samples of the continuous costates -- the covector mapping, made a
%   reusable library function. The principle: stationarity of the NLP
%   Lagrangian with respect to the states IS a discretization of the adjoint
%   equation lambda_dot = -dH/dx, so the defect multipliers sample the
%   costate. Which samples, with what scaling and sign, depends on the
%   transcription and the solver's conventions -- exactly what this function
%   owns, together with the verification that makes the output trustworthy.
%
%   WHAT IT RETURNS IS A SEED. Mapped costates are accurate to about 1e-3
%   regardless of care; single shooting amplifies that by ~1e3 over
%   multi-revolution arcs (measured: 12,700-558,000 km terminal misses).
%   Feed the output to multiple shooting (ms_tfmin); never treat it as a
%   shooting-quality solution.
%
%  ASSUMPTIONS / NOTES:
%
% • Station association is SCHEME-SPECIFIC: interval midpoints for
%   Hermite-Simpson, nodes for trapezoidal. Getting this wrong injects a
%   half-step phase error into every seed (found by external review).
% • CasADi callers: harvest multipliers from opti.lam_g indexed by the
%   defect rows. opti.dual returns duals in a CANONICALIZED orientation and
%   once caused a 10-60 degree primer misalignment.
% • The global sign of the multipliers is convention-dependent; it is fixed
%   empirically here by requiring the implied primer direction
%   -lambda_v/|lambda_v| to align with the solved thrust direction (majority
%   vote over the first stations). Problems without a thrust-like control
%   can skip the vote and resolve sign by their own physics.
% • On intervals where a PATH CONSTRAINT is active, inequality multipliers
%   contaminate the mapping: the result is the costate of the CONSTRAINED
%   problem. Pass activeFlags so the diagnostics can say so.
% • Lifted-final-time transcriptions: pass the row of multipliers paired
%   with t_f as lamTf; after sign resolution its median must equal +1 (it
%   is -H with the H(tf)=0, +1-normalized convention). Junk solutions
%   violate this by orders of magnitude -- it is the cheapest quality gate
%   this stage has.
%
%% Inputs:
%
%  spec                     struct
%   .scheme                 char                    'hermite-simpson' |
%                                                   'trapezoid'
%   .mu                     [ns x M]                Defect multipliers, one
%                                                   column per interval
%                                                   (M = N intervals for
%                                                   both schemes here)
%   .tNodes                 [1 x N+1]               Node times (physical);
%                                                   pass [] with .tf set for
%                                                   a uniform mesh
%   .tf                     double                  Final time (required if
%                                                   tNodes is empty)
%   .defectScale            char                    'unit' (defects written
%                                                   so mu ~ lambda directly;
%                                                   default) | 'h' (defects
%                                                   h-scaled: mu divided by
%                                                   the interval width)
%   .uDir                   [3 x M] (optional)      Solved thrust directions
%                                                   at the stations; enables
%                                                   the primer sign vote
%   .velRows                [1 x 3]                 Rows of mu holding
%                                                   lambda_v [default 4:6]
%   .nVote                  double                  Stations in the sign
%                                                   vote [default 10]
%   .lamTf                  [1 x M] (optional)      Multiplier row paired
%                                                   with lifted t_f; enables
%                                                   the lambda_t = +1 check
%   .activeFlags            [1 x M] (optional)      true where a path
%                                                   constraint is active on
%                                                   the interval
%
%% Outputs:
%
%  lam                      [ns x M]                Sign-resolved costate
%                                                   samples at the stations
%
%  tStations                [1 x M]                 Station times (interval
%                                                   midpoints for
%                                                   Hermite-Simpson, left
%                                                   nodes for trapezoid)
%
%  diag_                    struct                  .sign (+1/-1),
%                                                   .voteMargin (0..1; NaN
%                                                   if no vote), .lamT
%                                                   (median sign-resolved
%                                                   lambda_t; NaN if not
%                                                   given), .lamTOK
%                                                   (|lamT - 1| < 1e-3),
%                                                   .activeFraction,
%                                                   .notes {cellstr}
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: recover a stored campaign cell's costates, with a deliberately
   %flipped sign, and show the vote putting it back:
     RD = ['/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/', ...
           'direct/results'];
     CC = load(fullfile(RD,'dsweep_12x12_cells.mat'));
     c  = CC.CELLS{2,5};                     % the certified global-min cell
     s  = struct('scheme','hermite-simpson', 'mu', -c.lamDef(1:7,:), ...
                 'tNodes', c.tNodes, 'uDir', c.Um(1:3,1:size(c.lamDef,2)), ...
                 'lamTf', -c.lamDef(8,:));
     [lam, tS, dg] = duals_to_costates(s);
     fprintf('vote sign %+d (margin %.2f), lambda_t = %.6f (ok=%d)\n', ...
             dg.sign, dg.voteMargin, dg.lamT, dg.lamTOK);
     fprintf('lambda(0) approx [%s]\n', sprintf('%.3f ', lam(:,1)));
     return
end

d = @(f,v) fieldd(spec, f, v);
scheme  = lower(d('scheme', 'hermite-simpson'));
mu      = spec.mu;
[ns, M] = size(mu);
velRows = d('velRows', 4:6);
nVote   = d('nVote', 10);
notes   = {};

%% Station times (scheme-specific association):
tN = d('tNodes', []);
if isempty(tN)
    assert(isfield(spec,'tf') && ~isempty(spec.tf), ...
           'duals_to_costates:tf', 'tNodes empty: spec.tf is required');
    tN = linspace(0, spec.tf, M+1);
end
tN = tN(:).';
assert(numel(tN) == M+1, 'duals_to_costates:size', ...
       'tNodes (%d) must be one longer than mu columns (%d)', numel(tN), M);
switch scheme
    case 'hermite-simpson'
        tStations = tN(1:end-1) + diff(tN)/2;      % interval MIDPOINTS
    case 'trapezoid'
        tStations = tN(1:end-1);                   % nodes
        notes{end+1} = 'trapezoid: last-node costate not sampled';
    otherwise
        error('duals_to_costates:scheme', ...
              'scheme ''%s'' not implemented (hermite-simpson | trapezoid)', ...
              scheme);
end

%% Defect normalization:
switch lower(d('defectScale', 'unit'))
    case 'unit'
        % defects written so multipliers approximate lambda directly
    case 'h'
        mu = mu ./ diff(tN);                       % undo h-scaling
        notes{end+1} = 'h-scaled defects: multipliers divided by interval width';
    otherwise
        error('duals_to_costates:defectScale', 'unknown defectScale');
end

%% Global sign, by primer-vs-control vote where a control is supplied:
sg = 1;  voteMargin = NaN;
uDir = d('uDir', []);
if ~isempty(uDir)
    nv = min([nVote, M, size(uDir,2)]);
    vote = 0;
    for k = 1:nv
        pv = -mu(velRows,k) / max(norm(mu(velRows,k)), eps);
        ud =  uDir(:,k)     / max(norm(uDir(:,k)),     eps);
        vote = vote + sign(pv.' * ud);
    end
    if vote < 0, sg = -1; end
    voteMargin = abs(vote) / nv;
    if voteMargin < 0.6
        notes{end+1} = sprintf(['weak sign vote (margin %.2f): treat the ', ...
                                'mapping as suspect'], voteMargin);
    end
else
    notes{end+1} = 'no uDir supplied: sign NOT resolved (returned as-is)';
end
lam = sg * mu;

%% lambda_t consistency (lifted-final-time transcriptions):
lamT = NaN;  lamTOK = false;
lamTf = d('lamTf', []);
if ~isempty(lamTf)
    lamT = median(sg * lamTf);
    lamTOK = abs(lamT - 1) < 1e-3;
    if ~lamTOK
        notes{end+1} = sprintf(['lambda_t = %.4g (expected +1): solution or ', ...
                                'conventions are suspect'], lamT);
    end
end

%% Active-path-constraint contamination:
af = d('activeFlags', []);
activeFraction = 0;
if ~isempty(af)
    activeFraction = nnz(af) / numel(af);
    if activeFraction > 0
        notes{end+1} = sprintf(['%.0f%% of intervals have an active path ', ...
            'constraint: mapped costates there belong to the CONSTRAINED ', ...
            'problem'], 100*activeFraction);
    end
end

diag_ = struct('sign', sg, 'voteMargin', voteMargin, 'lamT', lamT, ...
               'lamTOK', lamTOK, 'activeFraction', activeFraction, ...
               'notes', {notes});
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
