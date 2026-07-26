function [sigma, X0, U0, tauf0, fp] = chain_rung_seed_tulip(src, tfNew, pNew, extraFp)
% CHAIN_RUNG_SEED_TULIP  Alias-free cross-rung warm start for the sundman engine.
%
% Takes a previous rung's CONVERGED Sundman solution, rescales ONLY the time
% row to tfNew (same-mesh scalar rescale -- controls stay attached to their own
% spatial nodes; no sigma- or t-interpolation anywhere, so the phase-aliasing
% channel never opens), and re-maps through the house no-resample
% sundman_seed_map: fresh tauf0 for the new rung (never reuse the source's),
% endpoints pinned exactly. (2026-07-21 triage C5/C6 ladder-prep.)
%
% INPUTS:
%   src     - loaded cache struct with .out (solver struct, X 8xM with t=row 8,
%             U 4xM), .rv0 [1x6], .rvf [1x6], and optionally .fp  [struct]
%   tfNew   - target transfer time for the new rung [ND, scalar]
%   pNew    - cr3bp_lt_params struct for the NEW rung (thrust differs!) [struct]
%   extraFp - optional extra fingerprint fields (e.g. .note, .epsMin) [struct]
% OUTPUTS:
%   sigma [M'x1], X0 [8xM'], U0 [4xM'], tauf0 [scalar] - warm start for
%   casadi_minfuel_sundman at the new rung;  fp - fingerprint with
%   .chainedFrom provenance [struct]
% REFERENCES:
%   [1] sundman_seed_map.m (the no-resample map); [2] spec sec 4.
if nargin < 4 || isempty(extraFp), extraFp = struct(); end
srcThrust = NaN;
if isfield(src,'fp') && isfield(src.fp,'thrustN'), srcThrust = src.fp.thrustN; end
if isnan(srcThrust), srcThrust = 0.025; end          % legacy caches are nominal
assert(abs(srcThrust - pNew.thrustN) > 1e-12, 'chain_rung_seed_tulip:sameThrust', ...
    'source and target thrust are both %.4g N -- chaining to the same rung is a caller bug', srcThrust);
X = src.out.X;  U = src.out.U;
assert(size(X,1) >= 8, 'chain_rung_seed_tulip:badState', 'need 8-state X with t=row 8');
X(8,:) = X(8,:) * (tfNew / X(8,end));                % time row only; spatial mesh untouched
cfg = minfuel_config();
[sigma, X0, U0, tauf0] = sundman_seed_map(X(1:7,:), U, tfNew, X(8,:).', ...
                                          cfg.pSund, pNew.muStar, src.rv0, src.rvf);
extraFp.chainedFrom = sprintf('T=%.4gN', srcThrust);
extraFp.tf = tfNew;
fp = cr3bp_fingerprint(pNew, extraFp);
end
