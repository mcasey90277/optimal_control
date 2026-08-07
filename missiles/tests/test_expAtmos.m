function test_expAtmos()
%% Purpose:
%
%  Verify the exponential atmosphere reproduces sea-level density, decays by
%  exactly one e-fold per scale height, vectorizes, and returns pressure and
%  sound speed matching values computed independently of the implementation.
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

%% Sea level returns rho0 exactly:
    [rho,P,T,a] = coorbital.atmos.expAtmos(0);
    assert(abs(rho - c.rho0) < 1e-12,'sea-level density wrong');

%% One scale height is one e-fold:
       [rhoH,~,~,~] = coorbital.atmos.expAtmos(c.Hscale);
    assert(abs(rhoH - c.rho0/exp(1)) < 1e-12,'scale-height decay wrong');

%% Sea-level pressure and sound speed against values worked out by hand, so a
%% formula error cannot hide behind the implementation's own arithmetic:
    assert(abs(T - 250)      < 1e-9,  'isothermal temperature wrong');
    assert(abs(P - 87909.92) < 1e-2,  'sea-level pressure wrong: got %.2f Pa',P);
    assert(abs(a - 316.9677) < 1e-3,  'sea-level sound speed wrong: got %.4f m/s',a);

%% The ideal gas law then holds at an altitude the hand values do not cover:
    [rho2,P2,T2,~] = coorbital.atmos.expAtmos(40000);
    assert(abs(P2 - rho2*c.Rair*T2) < 1e-6,'P does not satisfy the ideal gas law');

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
