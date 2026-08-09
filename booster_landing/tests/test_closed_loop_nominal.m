% TEST_CLOSED_LOOP_NOMINAL  Zero-dispersion closed loop must reproduce the
% guidance: touchdown within 1 m and 0.1 m/s of the pad-at-rest target.
% Also: with a 50 m initial position offset, the tracker must still land
% inside the pad radius (the whole point of feedback).
%
% ADDED ASSERTIONS (task-7 fix report round 2, 2026-08-08; STRICTER, not
% loosened -- per that review, a vertical arrest just above the pad can
% pass the miss/vtd checks above with an optimistically small vtd (vz=0
% is exactly the arrest condition) while never actually touching down.
% out0/out1.td.landed distinguishes a real z=0 touchdown from that
% arrest; both cases must genuinely land, not just clear the miss/vtd
% numbers.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
ctrl = tvlqr_design(sol, P);

out0 = sim_closed_loop(sol, ctrl, P, struct());
assert(out0.td.miss < 1.0,  'nominal miss %.2f m', out0.td.miss);
assert(out0.td.vtd  < 0.1 + 1e-9, 'nominal touchdown speed %.3f', out0.td.vtd);
assert(out0.td.landed, 'nominal did not land (stop=%s, alt=%.3f m) -- an arrest, not a touchdown', ...
       out0.td.stop, out0.td.alt);

dsp  = struct('dr0', [50; -30; 0]);
out1 = sim_closed_loop(sol, ctrl, P, dsp);
assert(out1.td.miss < P.pad_radius, 'dispersed miss %.1f m', out1.td.miss);
assert(out1.td.vtd  < P.vtd_max,    'dispersed vtd %.2f m/s', out1.td.vtd);
assert(out1.td.landed, 'dispersed did not land (stop=%s, alt=%.3f m) -- an arrest, not a touchdown', ...
       out1.td.stop, out1.td.alt);
fprintf('test_closed_loop_nominal PASS  (nom miss %.3f m, disp miss %.2f m)\n', ...
        out0.td.miss, out1.td.miss);
