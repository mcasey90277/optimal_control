function cert = table3_certified(thrustN)
% TABLE3_CERTIFIED  Certified Table-3 reproducer numbers for one thrust rung.
%
% Pure lookup: no file I/O, no solves, no side effects. Source of truth
% for the Table-3 reproducer engine's verify step (verify_row.m) -- these
% are the campaign's own certified numbers (see process/CAMPAIGN.md), not the
% published-paper numbers (those live in gergaud_row.m's revs_paper).
%
% INPUTS:
%   thrustN - max thrust level [N]; one of 10, 5, 2.5, 1, 0.5, 0.2, 0.1   [scalar]
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
%   [3] earth_elliptic_to_geo/process/DEEP_THRUST_LESSONS.md (0.2/0.1 N rungs,
%       certified 2026-07-20 via drivers/reproduce_deep_rung.m; numbers pulled
%       from the certified results/MEE_M2_0p{2,1}N.mat caches).
%   [4] earth_elliptic_to_geo/process/P0_SWITCH_MESH_CONVERGENCE.md.
%
% DEEP-RUNG SWITCH-COUNT CAVEAT (P0, 2026-07-21): the .switches values for the
% deep rungs (0.2 N = 823, 0.1 N = 1644) are 8-node/rev point estimates and are
% mesh-UNDER-resolved LOWER BOUNDS. Ref [4] refines 0.2 N to a converged band
% ~866+/-5 (823 is a ~5% undercount); mass and revs ARE mesh-converged. These
% fields are retained as the as-certified-at-8/rev values (what verify_row and
% the reproducer compare against); read deep-rung switch counts as bands, not
% exact integers.

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
    0.2, struct('m_f_kg',1377.2867,'switches',823, 'revs',346.726,'tfmin',1115.70, 'anchorSource','R0law'); ...
    0.1, struct('m_f_kg',1377.2882,'switches',1644,'revs',693.601,'tfmin',2231.40, 'anchorSource','R0law'); ...
};

for k = 1:size(rungs,1)
    if abs(rungs{k,1} - thrustN) < tol
        cert = rungs{k,2};
        cert.thrustN = thrustN;
        return;
    end
end

error('table3_certified:unknownThrust', ...
    'No certified Table-3 row for thrustN=%g (known: 10, 5, 2.5, 1, 0.5, 0.2, 0.1)', thrustN);

end
