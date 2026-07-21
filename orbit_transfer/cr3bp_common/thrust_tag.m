function tag = thrust_tag(thrustN)
% THRUST_TAG  Artifact filename token for a thrust rung.
%
% '' at the nominal 25 mN (all existing cache names stay byte-identical);
% otherwise '_T<mN>mN' with '.' -> 'p' (0.020 -> '_T20mN', 0.0325 -> '_T32p5mN').
%
% INPUTS:  thrustN - thrust [N, scalar]
% OUTPUTS: tag     - filename token [char]
% REFERENCES: [1] spec 2026-07-21-ladder-prep-design.md sec 2.
if abs(thrustN - 0.025) < 1e-12
    tag = '';
else
    mN = thrustN * 1000;
    s  = strrep(sprintf('%g', mN), '.', 'p');
    tag = sprintf('_T%smN', s);
end
end
