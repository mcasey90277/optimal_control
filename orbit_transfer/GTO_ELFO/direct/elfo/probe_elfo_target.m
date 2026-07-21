% PROBE_ELFO_TARGET  Build the ELFO in the ND barycentric rotating frame, pick an
% apolune rendezvous state, and report terminal stiffness vs the tulip target.
%
% The direct (PSR) pipeline is target-agnostic: minfuel_at_tf threads rv0/rvf
% from the seed file. So retargeting to an ELFO is entirely a matter of the
% rendezvous state rvf. This probe constructs the ELFO exactly as proj7 does
% (im_elfo_states.m: fromOrb about the Moon -> prop under CR3BP), picks apolune
% (highest, slowest, over the south pole -- the natural low-thrust capture
% point), and compares Earth/Moon distances against the tulip rvf and the GTO
% rv0 so we can gauge how much stiffer the terminal (lunar-perilune) arcs are.
%
% REFERENCES:
%   proj7/pipelines/gdop/im_elfo_states.m + im_elfo_optimum.m (ELFO construction)
%   gto_tulip_endpoints.m (tulip rvf, rv0)

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();

% --- CR3BP constants (match the tulip side exactly) -------------------------
mu    = 0.012150585609624;
Mmass = 5.9736e24 + 7.35e22;
lStar = 389703.264829278;
tStar = 382981.289129055;

% --- ELFO shape (im_elfo_optimum): sma,ecc,inc,argp; raan/M0 pick one member -
oe   = [12000, 0.69, 56.5, 90];      % [sma_km ecc inc_deg argp_deg]
raan = 0;                            % any member of the family is representative
M0   = 0;
nu0  = pumpkyn.util.mean2True(M0*pi/180, oe(2));
oev  = [oe(1), oe(2), oe(3)*pi/180, oe(4)*pi/180, raan*pi/180, nu0];
x0   = pumpkyn.cr3bp.fromOrb(0, oev, Mmass, mu, tStar, lStar, 2);   % ND baryc rot

% --- propagate one+ ELFO period and locate apolune / perilune ---------------
% ELFO period ~ 2*pi*sqrt(a^3/muMoon): a=12000 km, so ~33 h ~ 0.31 ND. Sample
% a bit past one period to bracket the first apolune cleanly.
tau  = linspace(0, 0.5, 6000)';
[~, xs] = pumpkyn.cr3bp.prop(tau, x0, mu);
rMoon = [1-mu, 0, 0];
dMoon = sqrt(sum((xs(:,1:3) - rMoon).^2, 2));
[dApo, iApo] = max(dMoon);
[dPer, iPer] = min(dMoon);
rvf_elfo = xs(iApo, 1:6);

% --- tulip target + GTO departure for comparison ----------------------------
p = cr3bp_lt_params(0.025, 15, 2100);
[rv0, rvf_tul] = gto_tulip_endpoints(p);

dfun = @(rv, c) norm(rv(1:3) - c) * lStar;    % ND -> km distance to center c
rE = [-mu 0 0];  rM = [1-mu 0 0];

fprintf('\n================= ELFO TARGET PROBE =================\n');
fprintf('ELFO shape: sma=%g km ecc=%.2f inc=%.1f deg argp=%.1f deg (raan=%g,M0=%g)\n',...
        oe(1),oe(2),oe(3),oe(4),raan,M0);
fprintf('ELFO period sampled over tau in [0,0.5] ND (%.1f h)\n', 0.5*tStar/3600);
fprintf('  apolune : %.1f km alt  (dist Moon %.1f km) at tau=%.4f\n', ...
        dApo*lStar-1737.4, dApo*lStar, tau(iApo));
fprintf('  perilune: %.1f km alt  (dist Moon %.1f km) at tau=%.4f\n', ...
        dPer*lStar-1737.4, dPer*lStar, tau(iPer));
fprintf('\n--- rendezvous state comparison (distances in km) ---\n');
fprintf('               |  dist Earth   dist Moon    speed(ND)\n');
fprintf('  GTO rv0      |  %9.1f   %9.1f    %.4f\n', dfun(rv0,rE),  dfun(rv0,rM),  norm(rv0(4:6)));
fprintf('  tulip rvf    |  %9.1f   %9.1f    %.4f\n', dfun(rvf_tul,rE), dfun(rvf_tul,rM), norm(rvf_tul(4:6)));
fprintf('  ELFO apolune |  %9.1f   %9.1f    %.4f\n', dfun(rvf_elfo,rE), dfun(rvf_elfo,rM), norm(rvf_elfo(4:6)));
fprintf('\nrvf_elfo (ND) = [% .8f % .8f % .8f % .8f % .8f % .8f]\n', rvf_elfo);
fprintf('rvf_tulip(ND) = [% .8f % .8f % .8f % .8f % .8f % .8f]\n', rvf_tul);
fprintf('||rvf_elfo - rvf_tulip|| = %.4f (ND state-space homotopy length)\n', norm(rvf_elfo-rvf_tul));
fprintf('====================================================\n\n');

save(fullfile(here,'results','elfo_target_probe.mat'), ...
     'rvf_elfo','rvf_tul','rv0','oe','raan','M0','dApo','dPer','lStar','mu');
