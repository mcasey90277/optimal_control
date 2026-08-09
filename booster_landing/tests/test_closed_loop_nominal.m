% TEST_CLOSED_LOOP_NOMINAL  Zero-dispersion closed loop must reproduce the
% guidance: a REAL touchdown (td.landed=true, not an arrest) within 1 m of
% the pad and within 0.1 m/s of the P.vf=-1.5 m/s descent target. Also:
% with a 50 m initial position offset, the tracker must still land inside
% the pad radius at <= P.vtd_max, and also genuinely land (not arrest).
%
% ADJUDICATED 2026-08-08 (task-7 fix report round 3): guidance terminal
% velocity changed from v(tf)=0 (singular under P.Tmin>weight -- see
% booster_params.m) to v(tf)=P.vf=[0;0;-1.5]. The |vtd-1.5| gate below
% (not vtd<0.1 outright) is the tracking-error check appropriate to that
% nonzero nominal descent speed -- the mission gate vtd<P.vtd_max=2.0 for
% the dispersed case is unchanged.
%
% ADDED ASSERTIONS (task-7 fix report round 2, 2026-08-08; STRICTER, not
% loosened -- per that review, a vertical arrest just above the pad can
% pass the miss/vtd checks above with an optimistically small vtd (vz=0
% is exactly the arrest condition) while never actually touching down.
% out0/out1.td.landed distinguishes a real z=0 touchdown from that
% arrest; both cases must genuinely land, not just clear the miss/vtd
% numbers. (Round 3's P.vf change removes the singular v=0-at-z=0
% endpoint that caused the arrest in the first place -- see
% booster_params.m -- so both blocks below are now expected to pass.)
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
ctrl = tvlqr_design(sol, P);

out0 = sim_closed_loop(sol, ctrl, P, struct());
assert(out0.td.miss < 1.0,  'nominal miss %.2f m', out0.td.miss);
assert(abs(out0.td.vtd - 1.5) < 0.1, 'nominal touchdown speed %.3f (target %.1f)', ...
       out0.td.vtd, -P.vf(3));
assert(out0.td.landed, 'nominal did not land (stop=%s, alt=%.3f m) -- an arrest, not a touchdown', ...
       out0.td.stop, out0.td.alt);

dsp  = struct('dr0', [50; -30; 0]);
out1 = sim_closed_loop(sol, ctrl, P, dsp);
assert(out1.td.miss < P.pad_radius, 'dispersed miss %.1f m', out1.td.miss);
assert(out1.td.vtd  < P.vtd_max,    'dispersed vtd %.2f m/s', out1.td.vtd);
assert(out1.td.landed, 'dispersed did not land (stop=%s, alt=%.3f m) -- an arrest, not a touchdown', ...
       out1.td.stop, out1.td.alt);
fprintf('test_closed_loop_nominal PASS  (nom miss %.3f m vtd %.3f m/s, disp miss %.2f m vtd %.3f m/s)\n', ...
        out0.td.miss, out0.td.vtd, out1.td.miss, out1.td.vtd);
