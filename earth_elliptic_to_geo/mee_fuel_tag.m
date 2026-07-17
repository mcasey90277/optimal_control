function tag = mee_fuel_tag(thrustN)
% MEE_FUEL_TAG  Canonical results/ cache tag for a thrust level's certified
% fuel (eps=0 homotopy) solution, e.g. run_transfer_mee.m's output. Single
% source of truth shared by run_mintime_mee.m (Stage A fuel-anchor lookup)
% and run_ladder.m (per-rung fuel-artifact reuse/production), so the two
% drivers can never disagree about where a given thrust's fuel result lives.
%
% Convention: integer thrusts -> 'MEE_M2_<N>N' (matches run_transfer_mee.m's
% own default cfg.tag='MEE_M2_10N' at thrustN=10, the Task-4 cross-
% formulation gate); non-integer thrusts -> decimal point replaced with 'p',
% e.g. thrustN=2.5 -> 'MEE_M2_2p5N'.
%
% INPUTS:  thrustN - max thrust [N, scalar]
% OUTPUTS: tag - cache tag stem [char], the file is resDir/[tag '.mat']
%
% REFERENCES: [1] earth_elliptic_to_geo/run_transfer_mee.m (producer).
%             [2] .superpowers/sdd/task-6-brief.md (Task 6, first consumer).
if abs(thrustN - round(thrustN)) < 1e-9
    tag = sprintf('MEE_M2_%dN', round(thrustN));
else
    s = strrep(sprintf('%g', thrustN), '.', 'p');
    tag = sprintf('MEE_M2_%sN', s);
end
end
