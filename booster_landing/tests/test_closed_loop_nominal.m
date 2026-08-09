% TEST_CLOSED_LOOP_NOMINAL  Zero-dispersion closed loop must reproduce the
% guidance: a REAL touchdown (td.landed=true, not an arrest) within 1 m of
% the pad, at or under the P.vtd_max mission gate. Also: with a 50 m
% initial position offset, the tracker must still land inside the pad
% radius at <= P.vtd_max, and also genuinely land (not arrest).
%
% Uses the PRODUCTION grid (N=P.N=60), not a coarse N=30 shortcut --
% ADAPTATION (task-7 fix report round 4, 2026-08-08): N=30 is exactly what
% caused two rounds of chasing a phantom "structural" tension. At N=30,
% hs_quad_ctrl's (pre-round-4) reconstruction dipped up to 18% below Tmin
% over 11.7% of the flight BETWEEN nodes/midpoints (invisible to every
% existing gate); the closed-loop sim's saturating clamp turned that into
% a real, grid-resolution-dependent open-loop altitude bias that TVLQR
% weight tuning was, in hindsight, fighting rather than a genuine
% closed-loop/physics limitation. hs_quad_ctrl is now annulus-feasible BY
% CONSTRUCTION at any grid (direction/magnitude reconstructed separately,
% magnitude clamped -- see that function's own note, and certify_pdg's new
% G2ff gate, which would have caught this from round 1), but the test
% still uses the production grid going forward since that is what the
% shipped defaults are tuned against.
%
% ADJUDICATED 2026-08-08 (task-7 fix report round 3): guidance terminal
% velocity changed from v(tf)=0 (singular under P.Tmin>weight -- see
% booster_params.m) to v(tf)=P.vf=[0;0;-1.5].
%
% GATE FIX (task-7 fix report round 4): the nominal block previously
% asserted |vtd-1.5|<0.1 -- a reviewer pointed out this is a hidden
% ~6.8 mm ALTITUDE requirement in disguise (dvz/dz ~ 14.7/s on the
% near-vertical terminal arc, so a 0.1 m/s vtd window maps to a
% millimeter-scale touchdown-altitude window that has nothing to do with
% the actual mission requirement). Replaced with the real mission gates:
% landed, miss<1m, vtd<=P.vtd_max -- vtd and terminal altitude are
% PRINTED as diagnostics (expected vtd ~1.5-1.9 m/s, comfortably under
% the 2.0 m/s gate) rather than asserted to a tight window.
%
% ADDED ASSERTIONS (task-7 fix report round 2, 2026-08-08; STRICTER, not
% loosened -- per that review, a vertical arrest just above the pad can
% pass the miss/vtd checks above with an optimistically small vtd (vz=0
% is exactly the arrest condition) while never actually touching down.
% out0/out1.td.landed distinguishes a real z=0 touchdown from that arrest.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', P.N));
ctrl = tvlqr_design(sol, P);

out0 = sim_closed_loop(sol, ctrl, P, struct());
fprintf('nominal diagnostics: vtd=%.3f m/s (target %.1f, gate <=%.1f)  alt=%.4f m  stop=%s\n', ...
        out0.td.vtd, -P.vf(3), P.vtd_max, out0.td.alt, out0.td.stop);
assert(out0.td.landed, 'nominal did not land (stop=%s, alt=%.3f m) -- an arrest, not a touchdown', ...
       out0.td.stop, out0.td.alt);
assert(out0.td.miss < 1.0,      'nominal miss %.2f m', out0.td.miss);
assert(out0.td.vtd  <= P.vtd_max, 'nominal touchdown speed %.3f m/s exceeds vtd_max %.1f', ...
       out0.td.vtd, P.vtd_max);

dsp  = struct('dr0', [50; -30; 0]);
out1 = sim_closed_loop(sol, ctrl, P, dsp);
fprintf('dispersed diagnostics: vtd=%.3f m/s  miss=%.3f m  alt=%.4f m  stop=%s\n', ...
        out1.td.vtd, out1.td.miss, out1.td.alt, out1.td.stop);
assert(out1.td.landed, 'dispersed did not land (stop=%s, alt=%.3f m) -- an arrest, not a touchdown', ...
       out1.td.stop, out1.td.alt);
assert(out1.td.miss < P.pad_radius, 'dispersed miss %.1f m', out1.td.miss);
assert(out1.td.vtd  < P.vtd_max,    'dispersed vtd %.2f m/s', out1.td.vtd);
fprintf('test_closed_loop_nominal PASS  (nom miss %.3f m vtd %.3f m/s, disp miss %.2f m vtd %.3f m/s)\n', ...
        out0.td.miss, out0.td.vtd, out1.td.miss, out1.td.vtd);
