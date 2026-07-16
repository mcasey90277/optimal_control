% TEST_SEED  Seed exactness: mass linearity, energy growth, small stencil defect,
% longitude-vs-throttle monotonicity, event stop at GEO energy.
p  = kepler_lt_params(10, 1500, 2000);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0, 0, pi, p.mu);   % coplanar
rv0 = [r0; v0];
% (a) fixed duration, full throttle
[sg, X0, U0, tauf0, inf1] = seed_2body(p, rv0, struct('sbar',1,'tDur',10,'N',400));
assert(abs(X0(7,end) - (1 - p.Tmax/p.c*10)) < 1e-8, 'mass not linear');
assert(all(diff(X0(8,:)) > 0) && abs(X0(8,end)-10) < 1e-9, 'time state bad');
assert(size(X0,1) == 9 && all(X0(9,:) == 1), 'cScale row missing');
% trapezoid defect of the mapped seed on the solver stencil must be small
dmax = seed_stencil_defect(sg, X0, U0, tauf0, p);
assert(dmax < 1e-2, 'seed defect too big: %.2e', dmax);
% (b) energy event stop: reaches GEO energy
[~, X1, ~, ~, inf2] = seed_2body(p, rv0, struct('sbar',1,'tDur',[],'N',200));
Eend = 0.5*norm(X1(4:6,end))^2 - 1/norm(X1(1:3,end));
assert(abs(Eend - (-0.5)) < 5e-3, 'did not stop at GEO energy');
% (c) winding decreases with throttle
[~,~,~,~, iLo] = seed_2body(p, rv0, struct('sbar',0.5,'tDur',30,'N',200));
[~,~,~,~, iHi] = seed_2body(p, rv0, struct('sbar',0.9,'tDur',30,'N',200));
assert(iLo.Larr > iHi.Larr, 'winding not monotone in sbar');
% (d) bisection hits a target longitude
Lt = inf2.Larr + 6;    % ask for ~1 rev more than the s=1 arrival
[~,~,~,~, iT] = seed_2body(p, rv0, struct('sbar',0.7,'tDur',1.5*inf2.tEnd,'N',200,'targetLf',Lt));
assert(abs(iT.Larr - Lt) < pi/2, 'targetLf bisection missed: %.2f vs %.2f', iT.Larr, Lt);
fprintf('test_seed: ALL PASS\n');

function dmax = seed_stencil_defect(sg, X0, U0, tauf0, p)
% Trapezoid defect of the seed under the solver's tau-stencil (numeric mirror).
N = numel(sg)-1;  dmax = 0;
f = zeros(9, N+1);
for k = 1:N+1
    rn  = norm(X0(1:3,k));
    fd  = lt2b_rhs_time(X0(1:8,k), U0(:,k), p);
    f(:,k) = [X0(9,k) * rn^p.pSund * fd; 0];
end
for k = 1:N
    d = X0(:,k+1) - X0(:,k) - (tauf0*(sg(k+1)-sg(k))/2)*(f(:,k)+f(:,k+1));
    dmax = max(dmax, max(abs(d)));
end
end
