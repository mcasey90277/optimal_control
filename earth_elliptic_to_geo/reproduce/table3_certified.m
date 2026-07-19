function cert = table3_certified(thrustN)
% TABLE3_CERTIFIED  Certified Table-3 reproducer numbers for one thrust rung.
%
% Pure lookup: no file I/O, no solves, no side effects. Source of truth
% for the Table-3 reproducer engine's verify step (verify_row.m) -- these
% are the campaign's own certified numbers (see process/CAMPAIGN.md), not the
% published-paper numbers (those live in gergaud_row.m's revs_paper).
%
% INPUTS:
%   thrustN - max thrust level [N]; one of 10, 5, 2.5, 1, 0.5             [scalar]
%
% OUTPUTS:
%   cert - struct with fields:
%     .thrustN     - echoes the input thrust level [N]                   [scalar]
%     .m_f_kg      - certified final mass [kg]                           [scalar]
%     .switches    - certified thrust-switch count                       [scalar]
%     .revs        - certified revolution count                          [scalar]
%     .tfmin       - certified min-time anchor, nondim time [TU]         [scalar]
%     .anchorSource - provenance of .tfmin: 'solved' (independently
%                    re-solved min-time anchor) or 'R0law' (the R0/T
%                    estimate; no anchor solve at this rung)              [char]
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, "Low Thrust Minimum-Fuel Orbital
%       Transfer: A Homotopic Approach," JGCD 27(6), 2004, Table 3.
%   [2] earth_elliptic_to_geo/process/CAMPAIGN.md (campaign record; certified
%       numbers harvested from here).

if nargin < 1
    error('table3_certified:badInput', 'thrustN is required');
end

tol = 1e-9;
rungs = { ...
    10,  struct('m_f_kg',1377.10, 'switches',19,  'revs',7.326,  'tfmin',22.2206, 'anchorSource','solved'); ...
    5,   struct('m_f_kg',1364.54, 'switches',32,  'revs',14.157, 'tfmin',44.6796, 'anchorSource','solved'); ...
    2.5, struct('m_f_kg',1369.79, 'switches',76,  'revs',27.841, 'tfmin',89.253,  'anchorSource','solved'); ...
    1,   struct('m_f_kg',1371.44, 'switches',171, 'revs',69.152, 'tfmin',223.808, 'anchorSource','solved'); ...
    0.5, struct('m_f_kg',1375.28, 'switches',362, 'revs',138.597,'tfmin',446.28,  'anchorSource','R0law'); ...
};

for k = 1:size(rungs,1)
    if abs(rungs{k,1} - thrustN) < tol
        cert = rungs{k,2};
        cert.thrustN = thrustN;
        return;
    end
end

error('table3_certified:unknownThrust', ...
    'No certified Table-3 row for thrustN=%g (known: 10, 5, 2.5, 1, 0.5)', thrustN);

end
