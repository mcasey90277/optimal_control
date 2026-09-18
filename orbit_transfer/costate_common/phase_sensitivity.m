function [dT_ds0, dT_dsf] = phase_sensitivity(lam0, lamF, dx0_ds, dxf_ds)
%% Purpose:
%
%   How much does the MINIMUM TIME of a transfer change when its departure
%   point, or its arrival point, slides along its orbit? The costates already
%   know. For a transfer from x_0(s_0) to x_f(s_f), each endpoint a point on
%   a curve parameterised by a phase s,
%
%       dT/ds_0 = + lam(0)   . dx_0/ds_0
%       dT/ds_f = - lam(t_f) . dx_f/ds_f
%
%   with lam the costates of the endpoint coordinates that are FIXED by the
%   boundary conditions (position and velocity; not the mass, whose initial
%   value does not depend on the phase and whose final value is free).
%
%   WHY. In the convention used throughout this library the Hamiltonian is
%   H = 1 + lam . f and the normal extremal has lam . f = -1, so lam(t) is the
%   gradient of the cost-to-go, lam = dV/dx. Moving the START by dx_0 changes
%   the time by lam(0) . dx_0 -- that is what "gradient of the cost-to-go"
%   means. At the END the target enters through the constraint
%   x(t_f) - x_f = 0 with multiplier nu = lam(t_f); moving the TARGET by dx_f
%   relaxes that constraint by -dx_f, hence the minus sign.
%
%   WHAT IT IS FOR. (1) A cross-check of a certified entry that no other gate
%   supplies: re-solve the transfer at s +/- delta, difference the two flight
%   times, and compare -- the re-solve never looks at a costate, so agreement
%   validates lam(0) and lam(t_f) independently (measured on the 70 mN
%   library: 4 to 5 digits). (2) The first-order condition of the FREE-phase
%   problem: a transfer is stationary with respect to a phase exactly where
%   its sensitivity vanishes. The catalog's entries are minima at FIXED
%   phases; this says which of them are candidates for the orbit-to-orbit
%   minimum, and which way the better neighbours lie.
%
%  ASSUMPTIONS / NOTES:
%
% • Minimum time, normal chart (lam . f = -1). For a cost J with running
%   cost L the same formulas give dJ/ds.
% • dx/ds is the derivative of the SAME endpoint map the solver used (the
%   phase closure), in the same phase units, so dT/ds comes out in time units
%   per unit phase. For a periodic orbit of period tau traversed by the flow,
%   dx/ds = tau * f_ballistic(x).
% • The test (tests/test_phase_sensitivity) uses the minimum-time single
%   integrator, T = |x_f - x_0|, where everything is known in closed form.
%
%% Inputs:
%
%  lam0                     [n x 1]                 costates of the fixed
%                                                   endpoint coordinates at t = 0
%  lamF                     [n x 1]                 the same at t = t_f
%  dx0_ds                   [n x 1]                 d(departure point)/d(phase)
%  dxf_ds                   [n x 1]                 d(arrival point)/d(phase)
%
%% Outputs:
%
%  dT_ds0                   double                  dT / d(departure phase)
%  dT_dsf                   double                  dT / d(arrival phase)
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%
%% References:
%   [1] Bryson & Ho, "Applied Optimal Control," 1975, Sec. 2.7-2.8 (the
%       costate as the sensitivity of the optimal cost to the state).
%% ------------------------ Begin Code Sequence ---------------------------

v = {lam0, lamF, dx0_ds, dxf_ds};
n = numel(lam0);
for k = 1:4
    assert(isnumeric(v{k}) && isvector(v{k}) && numel(v{k}) == n && all(isfinite(v{k})), ...
           'phase_sensitivity:input', 'the four inputs must be finite vectors of one length (%d)', n);
end
dT_ds0 =  sum(lam0(:).*dx0_ds(:));
dT_dsf = -sum(lamF(:).*dxf_ds(:));
end
