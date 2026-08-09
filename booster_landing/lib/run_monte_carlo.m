function mc = run_monte_carlo(sol, ctrl, P, opts)
% RUN_MONTE_CARLO  Dispersed closed-loop landing campaign.
%
% Draws initial-state, thrust-bias, Isp-bias (and, with drag on, wind)
% dispersions and runs sim_closed_loop per sample. Deterministic under
% P.seed (rng(P.seed) once, all draws via randn in a fixed order).
%
% TASK-7 EVOLUTION (binding, supersedes the task-8 brief): sim_closed_loop
% now classifies every run as touchdown/arrest/horizon (out.td.landed,
% out.td.stop). An arrest or horizon run is a FAILURE even when its
% miss/vtd/m look good -- an arrest's vtd is structurally optimistic (vz=0
% IS the arrest event, so vtd excludes the vertical-speed failure it
% actually represents). mc.ok therefore ALSO requires out.td.landed, and
% the failure-mode breakdown (n_landed/n_arrest/n_horizon) is tracked so a
% caller can tell a tracking failure from a genuine miss/vtd failure.
%
% INPUTS:
%   sol   - solve_pdg_colloc solution struct (Task 3)
%   ctrl  - tvlqr_design gain/feedforward struct (Task 6)
%   P     - booster_params struct
%   opts  - .Nrun [scalar, def 200]
%           .sig  [struct, 1-sigma dispersion magnitudes, any subset;
%                  defaults: .r0=[100;100;50] m, .v0=[10;10;10] m/s,
%                  .thrust=0.015 (1.5% of nominal), .isp=0.01 (1% of
%                  nominal), .wind=[10;10;0] m/s (drawn only when
%                  P.drag.on)]
%
% OUTPUTS:
%   mc - .land(Nrun x 2)   touchdown/terminal xy [m]
%        .vtd(Nrun x 1)    touchdown/terminal speed [m/s]
%        .mprop(Nrun x 1)  propellant remaining above P.mdry [kg]
%        .landed(Nrun x 1) logical, out.td.landed per run (real z=0
%                           touchdown vs. arrest/horizon -- see
%                           sim_closed_loop's ADAPTATION note)
%        .stop{Nrun x 1}   cellstr, out.td.stop per run ('touchdown' |
%                           'arrest' | 'horizon')
%        .ok(Nrun x 1)     logical: landed AND miss<P.pad_radius AND
%                           vtd<P.vtd_max AND m>=P.mdry
%        .success_rate     mean(mc.ok)
%        .n_landed/.n_arrest/.n_horizon  counts over Nrun (sum to Nrun)
%        .dr0/.dv0(Nrun x 3), .thrust_scale/.isp_scale(Nrun x 1) [,.wind
%          (Nrun x 3) when P.drag.on] -- the drawn dispersions themselves,
%          kept so a failure population's draw magnitudes are attributable
%          post hoc (e.g. "the failures are all >100 m lateral offsets").
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
%   [2] lib/sim_closed_loop.m -- .td.landed/.stop classification consumed
%       here; task-7 ADAPTATION note explains why an arrest cannot be
%       scored a success from .miss/.vtd alone.
if nargin < 4, opts = struct(); end
if ~isfield(opts,'Nrun'), opts.Nrun = 200; end
sig = struct('r0',[100;100;50], 'v0',[10;10;10], 'thrust',0.015, ...
             'isp',0.01, 'wind',[10;10;0]);
if isfield(opts,'sig')
    fn = fieldnames(opts.sig);
    for k = 1:numel(fn), sig.(fn{k}) = opts.sig.(fn{k}); end
end

rng(P.seed);
Nr = opts.Nrun;
mc.land         = zeros(Nr,2);
mc.vtd          = zeros(Nr,1);
mc.mprop        = zeros(Nr,1);
mc.landed       = false(Nr,1);
mc.stop         = cell(Nr,1);
mc.ok           = false(Nr,1);
mc.dr0          = zeros(Nr,3);
mc.dv0          = zeros(Nr,3);
mc.thrust_scale = zeros(Nr,1);
mc.isp_scale    = zeros(Nr,1);
if P.drag.on, mc.wind = zeros(Nr,3); end

for krun = 1:Nr
    d = struct('dr0', sig.r0 .* randn(3,1), ...
               'dv0', sig.v0 .* randn(3,1), ...
               'thrust_scale', 1 + sig.thrust*randn(), ...
               'isp_scale',    1 + sig.isp*randn());
    if P.drag.on, d.wind = sig.wind .* randn(3,1); end

    out = sim_closed_loop(sol, ctrl, P, d);

    mc.land(krun,:)        = out.td.r(1:2).';
    mc.vtd(krun)           = out.td.vtd;
    mc.mprop(krun)         = out.td.m - P.mdry;
    mc.landed(krun)        = out.td.landed;
    mc.stop{krun}          = out.td.stop;
    mc.ok(krun)            = out.td.landed && out.td.miss < P.pad_radius && ...
                              out.td.vtd < P.vtd_max && out.td.m >= P.mdry;
    mc.dr0(krun,:)         = d.dr0.';
    mc.dv0(krun,:)         = d.dv0.';
    mc.thrust_scale(krun)  = d.thrust_scale;
    mc.isp_scale(krun)     = d.isp_scale;
    if P.drag.on, mc.wind(krun,:) = d.wind.'; end

    if mod(krun, 25) == 0, fprintf('  MC %d/%d\n', krun, Nr); end
end

mc.success_rate = mean(mc.ok);
mc.n_landed     = nnz(mc.landed);
mc.n_arrest     = nnz(strcmp(mc.stop, 'arrest'));
mc.n_horizon    = nnz(strcmp(mc.stop, 'horizon'));
end
