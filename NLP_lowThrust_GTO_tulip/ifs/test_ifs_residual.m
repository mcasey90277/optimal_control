function test_ifs_residual()
% TEST_IFS_RESIDUAL  pack/unpack round-trip; residual squareness + zero on a
% self-consistent k=1 fixedState problem built by forward integration.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
p = cr3bp_lt_params(0.025, 15, 2100);  pSund = 1.5;
odeOpts = odeset('RelTol',1e-13,'AbsTol',1e-15);

% pack/unpack round trip
lam0 = (1:8).';  N = reshape(1:32, 16, 2);  tau = [3; 7];
Z = ifs_pack(lam0, N, tau);
assert(numel(Z) == 8 + 17*2, 'Z size 8+17k');
[l2, N2, t2] = ifs_unpack(Z, 2);
assert(isequal(l2,lam0) && isequal(N2,N) && isequal(t2,tau), 'round trip');

% Build a SELF-CONSISTENT k=1 fixedState problem by forward integration so the
% seed residual is ~0 (ground-truth-by-construction): pick a start state+costate,
% integrate a burn arc to tau1, then a coast arc to tauf, and adopt the endpoints
% as the node / terminal target. Choose costates so S changes sign at tau1.
%
% NOTE (deviation from the task-2 brief's literal Y0, disclosed in
% ifs-task-2-report.md): the brief's hand-picked costate tail gave
% S(N1)~=-12.5, not ~0, because it was never beta-scaled -- per the design doc
% (2026-07-11-ifs-design.md Sec.4), "the switching function's inhomogeneous
% '1' fixes the absolute costate scale...the dual-fit beta gives the right
% magnitude, not merely the direction." The costate tail below is the brief's
% direction ([-0.3 0.1 0.05 -0.7 0.2 0.03 -1.25 -0.4]) uniformly scaled by
% s=0.074005998841650, found via fzero so that S(N1)=0 to machine precision --
% i.e. tau1 is the TRUE switch time for this costate direction, matching the
% comment's stated intent.
Y0 = [0.9;0.0;-0.02; 0.05;1.3;-0.2; 1.0;0.0; ...
      -0.022201799652495; 0.007400599884165; 0.003700299942082; ...
      -0.051804199189155; 0.014801199768330; 0.002220179965249; ...
      -0.092507498552062; -0.029602399536660];
tau1 = 0.03;  tauf = 0.06;
[~,Ya] = ode113(@(s,y) ifs_eom(s,y,p.Tmax,p.c,p.muStar,pSund,1), [0 tau1], Y0, odeOpts);
N1 = Ya(end,:).';
[~,Yb] = ode113(@(s,y) ifs_eom(s,y,p.Tmax,p.c,p.muStar,pSund,0), [tau1 tauf], N1, odeOpts);
eEnd = Yb(end,:).';

prob = struct('rv0', Y0(1:6), 'm0', Y0(7), 't0', Y0(8), 'tau0', 0, ...
    'Tmax', p.Tmax, 'c', p.c, 'muStar', p.muStar, 'pSund', pSund, ...
    'tauf', tauf, 'k', 1, 'uArc', [1 0], 'termMode', 'fixedState', ...
    'termTarget', eEnd(1:8), 'odeOpts', odeOpts);
Zg = ifs_pack(Y0(9:16), N1, tau1);
R = ifs_residual(Zg, prob);
assert(numel(R) == 8 + 17*1, 'R square 8+17k');
assert(max(abs(R)) < 1e-8, 'ground-truth-by-construction residual ~0, got %.2e', max(abs(R)));

% rendezvous mode: residual has the right shape and is finite
probR = prob;  probR.termMode = 'rendezvous';
probR = rmfield(probR, 'termTarget');  probR.rvf = eEnd(1:6);  probR.tf = eEnd(8);
Rr = ifs_residual(Zg, probR);
assert(numel(Rr) == 8+17, 'rendezvous R square');
assert(all(isfinite(Rr)), 'rendezvous residual finite');
fprintf('ALL PASS (||R_fixedState||=%.2e)\n', max(abs(R)));
end
