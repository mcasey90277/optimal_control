% TEST_CLOSED_LOOP_NOMINAL  Zero-dispersion closed loop must reproduce the
% guidance: touchdown within 1 m and 0.1 m/s of the pad-at-rest target.
% Also: with a 50 m initial position offset, the tracker must still land
% inside the pad radius (the whole point of feedback).
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

dsp  = struct('dr0', [50; -30; 0]);
out1 = sim_closed_loop(sol, ctrl, P, dsp);
assert(out1.td.miss < P.pad_radius, 'dispersed miss %.1f m', out1.td.miss);
assert(out1.td.vtd  < P.vtd_max,    'dispersed vtd %.2f m/s', out1.td.vtd);
fprintf('test_closed_loop_nominal PASS  (nom miss %.3f m, disp miss %.2f m)\n', ...
        out0.td.miss, out1.td.miss);
