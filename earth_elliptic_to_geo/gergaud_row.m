function row = gergaud_row(inp)
% GERGAUD_ROW  Assemble one Table-3-style row (Haberkorn-Martinon-Gergaud,
% JGCD 27(6), 2004) for a min-fuel Earth-elliptic-to-GEO low-thrust result.
%
% Pure function: no file I/O, no solves, no side effects. Takes the raw
% fields of an already-solved/loaded result and derives propellant mass,
% delta-V (rocket equation), tfmin in hours, and the paper's own published
% Table-3 revolution count at the matching thrust level. All unit
% conversions go through kepler_lt_params (single source of physical
% constants) -- nothing here is a re-hardcoded magic number.
%
% INPUTS:
%   inp - struct with fields:
%     .thrustN   - max thrust [N] (paper cases: 10, 5, 2.5, 1, 0.5)        [scalar]
%     .tfmin_ND  - min-time anchor, nondim time [TU]                      [scalar]
%     .ctf       - tf/tfmin ratio used for this fuel solve                [scalar]
%     .tf_ND     - transfer time actually solved, nondim [TU]             [scalar]
%     .m_f_kg    - final mass [kg]                                        [scalar]
%     .switches  - number of thrust switches (bang-bang structure count)  [scalar]
%     .revs      - number of revolutions (ours)                          [scalar]
%     .edge      - control constraint edge/saturation fraction [0,1]      [scalar]
%     .incl_deg  - terminal inclination [deg]                             [scalar]
%     .defect    - collocation/shooting defect residual                   [scalar]
%     .certified - true if this result is independently certified        [logical]
%     .note      - free-text caveat/provenance note (e.g. "not attained") [char]
%     .m0kg      - OPTIONAL initial mass [kg] (default 1500)              [scalar]
%     .ispS      - OPTIONAL specific impulse [s] (default 2000)           [scalar]
%
% OUTPUTS:
%   row - struct: all of inp's fields passed through (plus resolved .m0kg
%         and .ispS), plus:
%     .prop_kg    - propellant consumed, m0kg - m_f_kg [kg]                    [scalar]
%     .dV_kms     - rocket-equation delta-V, c*log(m0kg/m_f_kg)*VU_kms [km/s]  [scalar]
%     .tfmin_h    - tfmin_ND converted to hours                               [scalar]
%     .revs_paper - paper's own published Table-3 revolution count at this
%                   thrust level (10->7.5, 5->15, 2.5->30, 1->74.5,
%                   0.5->149; NaN if thrustN matches none of the five)        [scalar]
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, "Low Thrust Minimum-Fuel Orbital
%       Transfer: A Homotopic Approach," JGCD 27(6), 2004, Table 3.
%   [2] kepler_lt_params.m (canonical unit conversions, same folder).

if nargin < 1 || ~isstruct(inp)
    error('gergaud_row:badInput', 'inp must be a struct');
end

if isfield(inp, 'm0kg') && ~isempty(inp.m0kg)
    m0kg = inp.m0kg;
else
    m0kg = 1500;
end
if isfield(inp, 'ispS') && ~isempty(inp.ispS)
    ispS = inp.ispS;
else
    ispS = 2000;
end

p = kepler_lt_params(inp.thrustN, m0kg, ispS);

row      = inp;            % pass through all input fields verbatim
row.m0kg = m0kg;
row.ispS = ispS;

row.prop_kg    = m0kg - inp.m_f_kg;
row.dV_kms     = p.c * log(m0kg / inp.m_f_kg) * p.VU_kms;
row.tfmin_h    = inp.tfmin_ND * p.TU_s / 3600;
row.revs_paper = lookup_paper_revs(inp.thrustN);

end

function rp = lookup_paper_revs(thrustN)
% Paper Table 3 (Haberkorn-Martinon-Gergaud 2004) published revolution
% counts, keyed by thrust level [N]. NaN for any thrustN off the ladder.
paperT    = [10, 5, 2.5, 1, 0.5];
paperRevs = [7.5, 15, 30, 74.5, 149];
tol = 1e-9;
idx = find(abs(paperT - thrustN) < tol, 1);
if isempty(idx)
    rp = NaN;
else
    rp = paperRevs(idx);
end
end
