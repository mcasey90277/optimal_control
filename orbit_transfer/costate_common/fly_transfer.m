function F = fly_transfer(z8, rv0, rvf, phys)
%% Purpose:
%
%   Fly converged min-time costates ONCE and hand back everything the
%   consumers then recompute: the flight itself, its admissibility, where it
%   actually arrived, and the mass and Delta-V that follow from it.
%
%   transfer_study section 5, certify_root, audit_phase_catalog and the
%   witness flight in verify_with_pumpkyn each spelled this out. ONE flight
%   object also answers the question a reader should never have to ask --
%   which trajectory owns a reported number.
%
%  ASSUMPTIONS / NOTES:
%
% • All-burn min-time flight (pumpkyn tfMinProp), augmented state
%   [r; v; m; lambda] with mass fraction 1 at t = 0.
%
% • It REPORTS admissibility, it does not throw on it: a script asserts, a
%   certifier returns a named reason, and this cannot know which it serves.
%
% • Delta-V is the rocket equation c*log(1/mf), the same one
%   `catalog_schema`'s `deltav_from_mf` derivation applies to catalog
%   entries. The schema's form needs a whole catalog struct, so the formula
%   is written here rather than faked through it; `tests/test_fly_transfer`
%   asserts the two agree on the real catalog.
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lambda0(7); tf]
%  rv0, rvf                 [6 x 1]                 Departure and arrival
%                                                   states (rvf is used for
%                                                   the flown miss only)
%  phys                     struct                  .Tnd .cnd .muStar
%                                                   .lStar .tStar .m0kg
%
%% Outputs:
%
%  F                        struct                  .t .Y (the flight)
%                                                   .z8 .rv0 .rvf .Tnd .cnd
%                                                   .muStar .lStar .tStar
%                                                   .m0kg .nSamples
%                                                   .admissibility
%                                                   (validate_flight)
%                                                   .mf .finalMassKg
%                                                   .propellantKg .dvKms
%                                                   .tfNd .tfDays
%                                                   .flyKm .flyVms
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if numel(z8) ~= 8
    error('fly_transfer:z8', 'z8 must be [lambda0(7); tf] (8 elements), got %d', numel(z8));
end
           need = {'Tnd', 'cnd', 'muStar', 'lStar', 'tStar', 'm0kg'};
for k = 1:numel(need)
    if ~isfield(phys, need{k}) || isempty(phys.(need{k}))
        error('fly_transfer:physics', 'phys.%s is required to fly a transfer', need{k});
    end
end
            z8 = z8(:);   rv0 = rv0(:);   rvf = rvf(:);

%% ONE flight:
        [t, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6); 1; z8(1:7)], ...
                                         phys.Tnd, phys.cnd, phys.muStar);
             F = struct('t', t, 'Y', Y, 'z8', z8, 'rv0', rv0(1:6), 'rvf', rvf(1:6), ...
                        'Tnd', phys.Tnd, 'cnd', phys.cnd, 'muStar', phys.muStar, ...
                        'lStar', phys.lStar, 'tStar', phys.tStar, 'm0kg', phys.m0kg, ...
                        'nSamples', numel(t));

%% Its admissibility, through the ONE validator the certifier also uses:
F.admissibility = validate_flight(t, Y, z8(8), phys.Tnd, phys.cnd, phys.muStar, phys.lStar);

%% What every consumer derives from it:
          F.mf = Y(end,7);
 F.finalMassKg = F.mf*phys.m0kg;
F.propellantKg = (1 - F.mf)*phys.m0kg;
       F.dvKms = phys.cnd*log(1/F.mf)*phys.lStar/phys.tStar;
        F.tfNd = z8(8);
      F.tfDays = z8(8)*phys.tStar/86400;
       F.flyKm = norm(Y(end,1:3).' - rvf(1:3))*phys.lStar;
      F.flyVms = norm(Y(end,4:6).' - rvf(4:6))*phys.lStar/phys.tStar*1000;
end
