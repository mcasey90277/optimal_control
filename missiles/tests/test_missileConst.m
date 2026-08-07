function test_missileConst()
%% Purpose:
%
%  Verify the constants struct exposes every field the library depends on,
%  with values in SI and of the right order of magnitude.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();

%% Every required field is present:
          required = {'muE','rE','omegaE','rho0','Hscale','T0','Rair','gamAir','g0'};
    for k = 1:numel(required)
        assert(isfield(c,required{k}), ...
            'missileConst is missing field %s',required{k});
    end

%% Values are physically sensible:
    assert(abs(c.muE  - 3.986004418e14) < 1e6,   'muE out of range');
    assert(abs(c.rE   - 6378137)        < 1,     'rE out of range');
    assert(abs(c.g0   - 9.80665)        < 1e-3,  'g0 out of range');
    assert(abs(c.rho0 - 1.225)          < 1e-3,  'rho0 out of range');
%% Hscale is PINNED, not range-checked. The Allen-Eggers reference
%% aMax = Ve^2|sin(gammaE)|/(2 e Hscale) divides by the same scale height the
%% atmosphere multiplies by, so that comparison is structurally blind to a wrong
%% Hscale -- both sides move together. This assertion is the only thing standing
%% between a mistyped scale height and a validation suite that still reports green:
    assert(abs(c.Hscale - 7200)         < 1,     'Hscale out of range');
    assert(abs(c.gamAir - 1.4)          < 1e-6,  'gamAir out of range');
    assert(abs(c.T0   - 250)            < 1e-9,  'T0 out of range');
    assert(abs(c.Rair - 287.053)        < 1e-3,  'Rair out of range');

%% Surface gravity from muE and rE agrees with g0 to better than 0.5 percent:
              gSurf = c.muE/c.rE^2;
    assert(abs(gSurf - c.g0)/c.g0 < 5e-3, ...
        'mu/r^2 = %.4f disagrees with g0 = %.4f',gSurf,c.g0);
end
