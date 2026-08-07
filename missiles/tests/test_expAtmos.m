function test_expAtmos()
%% Purpose:
%
%  Verify the exponential atmosphere reproduces sea-level density, decays by
%  exactly one e-fold per scale height, vectorizes, and returns a sound speed
%  consistent with its own temperature.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();

%% Sea level returns rho0 exactly:
    [rho,P,T,a] = coorbital.atmos.expAtmos(0);
    assert(abs(rho - c.rho0) < 1e-12,'sea-level density wrong');

%% One scale height is one e-fold:
       [rhoH,~,~,~] = coorbital.atmos.expAtmos(c.Hscale);
    assert(abs(rhoH - c.rho0/exp(1)) < 1e-12,'scale-height decay wrong');

%% Ideal gas law is self-consistent:
    assert(abs(P - rho*c.Rair*T) < 1e-6,'P does not satisfy the ideal gas law');

%% Sound speed matches its own temperature:
    assert(abs(a - sqrt(c.gamAir*c.Rair*T)) < 1e-9,'sound speed inconsistent');

%% Vectorizes over a column of altitudes:
              hVec = (0:10000:100000)';
    [rhoV,PV,TV,aV] = coorbital.atmos.expAtmos(hVec);
    assert(isequal(size(rhoV),size(hVec)),'rho not [N x 1]');
    assert(isequal(size(PV),  size(hVec)),'P not [N x 1]');
    assert(isequal(size(TV),  size(hVec)),'T not [N x 1]');
    assert(isequal(size(aV),  size(hVec)),'a not [N x 1]');
    assert(all(diff(rhoV) < 0),'density must decrease monotonically');

%% Negative altitude does not produce a non-finite value:
       [rhoNeg,~,~,~] = coorbital.atmos.expAtmos(-100);
    assert(isfinite(rhoNeg),'negative altitude gave a non-finite density');
end
