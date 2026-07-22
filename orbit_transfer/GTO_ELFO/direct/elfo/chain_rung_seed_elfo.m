function [seedS, fp] = chain_rung_seed_elfo(seedS, pNew, extraFp)
% CHAIN_RUNG_SEED_ELFO  Cross-rung warm start for the freetf engine (trivial by
% design: the cScale slack state decouples the clock, so a thrust rung reuses
% the source X,U,tauf0 on the SAME nodes -- no interpolation, no aliasing).
% Adds the new rung's fingerprint + chainedFrom provenance; refuses a
% same-thrust chain. (2026-07-21 ladder-prep, spec sec 4.)
%
% INPUTS:
%   seedS   - loaded seed struct (.X [9xM] or [8xM], .U [4xM], .tauf0, .sigma,
%             .rv0, .rvf, optional .fp) [struct]
%   pNew    - cr3bp_lt_params for the NEW rung [struct]
%   extraFp - optional extra fingerprint fields [struct]
% OUTPUTS:
%   seedS - unchanged seed struct (pass-through); fp - new-rung fingerprint
% REFERENCES: [1] casadi_energy_freetf.m (cScale mechanics); [2] spec sec 4.
if nargin < 3 || isempty(extraFp), extraFp = struct(); end
srcThrust = 0.025;                                    % legacy caches are nominal
if isfield(seedS,'fp') && isfield(seedS.fp,'thrustN'), srcThrust = seedS.fp.thrustN; end
assert(abs(srcThrust - pNew.thrustN) > 1e-12, 'chain_rung_seed_elfo:sameThrust', ...
    'source and target thrust are both %.4g N', srcThrust);
extraFp.chainedFrom = sprintf('T=%.4gN', srcThrust);
if isfield(seedS,'X'), extraFp.tf = seedS.X(8,end); end
fp = cr3bp_fingerprint(pNew, extraFp);
end
