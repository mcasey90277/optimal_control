### Task 9: Visualization (`plot_pdg_solution`, `plot_footprint`, `movie_landing`)

**Files:**
- Create: `viz/plot_pdg_solution.m`
- Create: `viz/plot_footprint.m`
- Create: `viz/movie_landing.m`
- Create: `tests/test_viz_smoke.m`

**Interfaces:**
- Consumes: `sol` (either solver), `mc`, `out` (closed-loop trace).
- Produces:
  - `fig = plot_pdg_solution(solC, solV, outfile)` — 2×2: (a) 3D trajectory + glideslope cone + thrust arrows every 5th node, both solutions overlaid; (b) throttle `‖T‖/Tmax` vs t (the max–min–max money plot), both solvers; (c) mass vs t; (d) speed vs t. `exportgraphics` to `outfile` PNG at 200 dpi if given.
  - `fig = plot_footprint(mc, P, outfile)` — landing scatter colored by success, pad-radius circle, 3σ ellipse (from `cov(mc.land)`), annotation with success rate and vtd stats.
  - `movie_landing(out, sol, P, outfile)` — house polished-graphics MP4: left panel booster altitude view (marker + thrust-vector arrow scaled by throttle, ground line, pad), right panel throttle trace with moving cursor; **frame size forced to 1280×720** (÷16 rule — the diagonal-streak H.264 lesson), 30 fps, `VideoWriter` MPEG-4. Loop over a fixed 12 s duration resampled from `out.t`.
- All three follow matlab-polished-graphics: fixed axis limits computed once (no autoscale jitter), title fixed-width, no emoji, colorblind-safe two-color scheme (`[0 0.447 0.741]` colloc / `[0.85 0.325 0.098]` convex).

- [ ] **Step 1: Write the failing smoke test**

`tests/test_viz_smoke.m`:

```matlab
% TEST_VIZ_SMOKE  Viz functions execute headless and write files.
% Not a beauty contest -- existence + nonzero size only. Movie smoke uses
% 2 seconds of frames to stay fast.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
solC = solve_pdg_colloc(P, struct('N', 30));
solV = solve_pdg_convex(P, struct('tf', solC.tf, 'Nconv', 60));
ctrl = tvlqr_design(solC, P);
out  = sim_closed_loop(solC, ctrl, P, struct());
mc   = run_monte_carlo(solC, ctrl, P, struct('Nrun', 8));

od = fullfile(tempdir, 'bl_viz_smoke');  if ~exist(od,'dir'), mkdir(od); end
plot_pdg_solution(solC, solV, fullfile(od, 'sol.png'));
plot_footprint(mc, P, fullfile(od, 'fp.png'));
movie_landing(out, solC, P, fullfile(od, 'mov.mp4'), struct('duration', 2));
for f = {'sol.png','fp.png','mov.mp4'}
    d = dir(fullfile(od, f{1}));
    assert(~isempty(d) && d.bytes > 1e3, 'missing/empty %s', f{1});
end
fprintf('test_viz_smoke PASS\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_viz_smoke"`
Expected: FAIL — functions undefined.

- [ ] **Step 3: Implement the three viz functions**

Follow the interface bullets above exactly; skeleton for the movie (the layout/streak-proofing parts that must not be improvised):

```matlab
function movie_landing(out, sol, P, outfile, opts)
% MOVIE_LANDING  Landing movie: booster + thrust vector | throttle trace.
% Frame size locked to 1280x720 (divisible by 16 -- H.264 shear guard).
% INPUTS: out (sim_closed_loop), sol (guidance), P, outfile .mp4,
%         opts.duration [s of playback, def 12], opts.fps [def 30]
if nargin < 5, opts = struct(); end
if ~isfield(opts,'duration'), opts.duration = 12; end
if ~isfield(opts,'fps'),      opts.fps = 30;      end

fig = figure('Visible','off','Position',[100 100 1280 720],'Color','w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
axT = nexttile(tl);  axR = nexttile(tl);

%% Fixed limits computed ONCE (no autoscale jitter):
xmax = max([abs(out.X(:,1)); abs(out.X(:,2)); 50]);
zmax = max(out.X(:,3)) * 1.05;
Tmag = sqrt(sum(out.Tcmd.^2, 2));

nf = round(opts.duration * opts.fps);
ts = interp1(linspace(0,1,numel(out.t)), out.t, linspace(0,1,nf));
vw = VideoWriter(outfile, 'MPEG-4');  vw.FrameRate = opts.fps;  open(vw);
for k = 1:nf
    xk = interp1(out.t, out.X, ts(k)).';
    Tk = interp1(out.t, out.Tcmd, ts(k)).';
    % left: altitude view (downrange = sqrt(x^2+y^2) signed by x), booster
    % marker, thrust arrow DOWN-scaled by throttle, pad + ground line
    cla(axT);  hold(axT, 'on');
    plot(axT, sqrt(sum(out.X(:,1:2).^2,2)).*sign(out.X(:,1)+eps), ...
         out.X(:,3), '-', 'Color', [0.7 0.7 0.7]);
    dk = sqrt(sum(xk(1:2).^2)) * sign(xk(1)+eps);
    arr = -0.15 * zmax * Tk / P.Tmax;    % plume points opposite thrust
    quiver(axT, dk, xk(3), arr(1), arr(3), 0, 'r', 'LineWidth', 2, ...
           'MaxHeadSize', 0.5);
    plot(axT, dk, xk(3), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 9);
    plot(axT, [-P.pad_radius P.pad_radius], [0 0], 'g-', 'LineWidth', 4);
    xlim(axT, [-xmax xmax]*1.1);  ylim(axT, [-0.02*zmax zmax]);
    title(axT, sprintf('t = %6.2f s   alt = %7.1f m', ts(k), xk(3)));
    xlabel(axT, 'downrange [m]');  ylabel(axT, 'altitude [m]');
    % right: throttle trace + cursor
    cla(axR);  hold(axR, 'on');
    plot(axR, out.t, Tmag/P.Tmax, 'b-');
    yline(axR, P.Tmin/P.Tmax, 'k--');  yline(axR, 1, 'k--');
    xline(axR, ts(k), 'r-');
    ylim(axR, [0 1.1]);  xlim(axR, [0 out.t(end)]);
    xlabel(axR, 'time [s]');  ylabel(axR, 'throttle T/Tmax');
    title(axR, 'throttle');
    frame = getframe(fig);
    frame.cdata = imresize(frame.cdata, [720 1280]);   % divide-by-16 lock
    writeVideo(vw, frame);
end
close(vw);  close(fig);
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_viz_smoke PASS`. Then open the smoke PNGs once (`open /tmp/.../sol.png` via the tempdir path printed) and eyeball: throttle plot must visibly show max–min–max.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/viz/plot_pdg_solution.m booster_landing/viz/plot_footprint.m \
        booster_landing/viz/movie_landing.m booster_landing/tests/test_viz_smoke.m
git commit -m "booster_landing: solution/footprint plots + landing movie (1280x720 lock)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

