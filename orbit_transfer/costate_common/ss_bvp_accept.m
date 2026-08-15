function [z, info] = ss_bvp_accept(prob, y1, tf, opts)
%% Purpose:
%
%   GENERIC SINGLE-SHOOTING ACCEPTANCE GATE -- the family- and cost-agnostic
%   form of the catalog pipeline's third gate ("feed the refined entry to
%   the single-shooting solver; PASS = it comes back unchanged in a few
%   iterations"). For min-time the independent solver is pumpkyn's tfMin;
%   for every other cost (min-energy, min-fuel) there is no pumpkyn twin, so
%   this harness plays the role on the same three-closure prob that ms_bvp
%   consumes. Mechanically it IS ms_bvp with K = 1 -- one full-arc
%   propagation, unknowns y1(freeIdx0) (+ tf when free), terminal residual
%   only -- which is exactly single shooting; the value added is the
%   acceptance semantics: the residual AT the seed before any step
%   (normR0: "is it already a root?"), the distance moved (dz), and the
%   verdict accepted = converged AND dz < tolDz.
%
%   TOLERANCE IS THE SINGLE-SHOOTING FLOOR, NOT ms_bvp's. A full-arc CR3BP
%   propagation at RelTol 1e-10 amplifies integrator round-off by the STM
%   growth (1e3-1e4 over 4-7 ND), so the single-shooting residual is only
%   defined to ~1e-7: measured on the min-energy pilot, the SAME lambda_0
%   propagated with and without the variational equations lands 6e-7 apart.
%   Hence the default .tolR here is 1e-6 (documented floor; pumpkyn tfMin's
%   own gate is likewise judged by |dz|, not by a 1e-10 residual). The tight
%   certification is the ms residual + the flown arrival; this gate certifies
%   USABILITY in a single-shooting workflow.
%
%  ASSUMPTIONS / NOTES:
%
% • Single shooting over the full arc is exactly the sensitivity-amplifying
%   solve ms_bvp exists to avoid -- which is the point: an entry that
%   survives it unchanged is a root of the single-shooting equations and
%   usable in any single-shooting workflow with zero adaptation.
% • Same prob contract as ms_bvp (prop must throw on collapse; terminal
%   supplies nf (+1 free-tf) conditions).
%
%% Inputs:
%
%  prob                     struct                  As ms_bvp (.ny,
%                                                   .freeIdx0, .prop, .rhs,
%                                                   .terminal)
%
%  y1                       [ny x 1]                Full initial state:
%                                                   fixed part + the
%                                                   refined free part
%
%  tf                       double                  Final time (fixed, or
%                                                   the refined value when
%                                                   free)
%
%  opts                     struct (optional)       .fixedTf [false],
%                                                   .tolDz [1e-6], plus any
%                                                   ms_bvp option (.tolR
%                                                   [1e-6 -- the single-
%                                                   shooting floor, see
%                                                   above], .maxIter [50],
%                                                   .wallSec [120], ...)
%
%% Outputs:
%
%  z                        [nf x 1] or [nf+1 x 1]  Returned unknowns
%                                                   [y1(freeIdx0)] (+ tf
%                                                   when free)
%
%  info                     struct                  .accepted, .dz (inf-norm
%                                                   move), .normR0 (residual
%                                                   at the seed, no-STM
%                                                   flight), .normR,
%                                                   .converged (at .tolR),
%                                                   .tolR, .iters, .wall
%
%% Revision History:
%  M. Casey                                                   (c) 08/14/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: oscillator, fixed tf = pi/2, seed exactly at the root y'(0) = 1.
     A_ = [0 1; -1 0];
     pr_ = struct('ny',2, 'freeIdx0',2, ...
         'prop', @(dt,y0,nS) demo_prop(A_,dt,y0,nS), 'rhs', @(y) A_*y, ...
         'terminal', @(y,nJ) deal(y(1)-1, [1 0]));
     [z_, inf_] = ss_bvp_accept(pr_, [0; 1], pi/2, struct('fixedTf', true));
     fprintf('demo: z = %.9f, normR0 = %.1e, dz = %.1e, accepted = %d\n', ...
             z_(1), inf_.normR0, inf_.dz, inf_.accepted);
     return
end

if nargin < 4, opts = struct(); end
fixedTf = isfield(opts, 'fixedTf') && opts.fixedTf;
tolDz   = 1e-6;
if isfield(opts, 'tolDz'), tolDz = opts.tolDz; end
mo = opts;
if isfield(mo, 'tolDz'), mo = rmfield(mo, 'tolDz'); end
if ~isfield(mo, 'maxIter'), mo.maxIter = 50; end
if ~isfield(mo, 'wallSec'), mo.wallSec = 120; end
if ~isfield(mo, 'tolR'),    mo.tolR    = 1e-6; end    % single-shooting floor
mo.fixedTf = fixedTf;

y1 = y1(:);
fi0 = prob.freeIdx0(:)';
zIn = y1(fi0);
if ~fixedTf, zIn(end+1) = tf; end

% Residual AT the seed: one propagation, no solve (K = 1 => the residual is
% just the terminal condition at the propagated endpoint).
yh = prob.prop(tf, y1, false);
[g0, ~] = prob.terminal(yh, false);      % two outputs: the ms_bvp contract
normR0 = max(abs(g0));

seed = struct('tf', tf, 'tGrid', [0 tf], 'Y', [y1, yh]);
[p, mi] = ms_bvp(prob, seed, mo);
z = p;
dz = max(abs(z - zIn));
info = struct('accepted', mi.converged && dz < tolDz, 'dz', dz, ...
              'normR0', normR0, 'normR', mi.normR, 'tolR', mo.tolR, ...
              'converged', mi.converged, 'iters', mi.iters, 'wall', mi.wall);
end

% ------------------------------------------------------------------------
function [yh, PHI] = demo_prop(A, dt, y0, needSTM)
% DEMO_PROP  Exact propagator for the linear demo.  INPUTS: A [2x2]; dt;
% y0 [2x1]; needSTM logical.  OUTPUTS: yh [2x1]; PHI [2x2] or [].
E = expm(A*dt);
yh = E*y0;
if needSTM, PHI = E; else, PHI = []; end
end
