function row = run_gergaud(opts)
% RUN_GERGAUD  Front-door entry point for one row of Gergaud-Haberkorn-
% Martinon (JGCD 27(6), 2004) Table 3 -- a min-fuel low-thrust transfer from
% an elliptic, inclined start orbit to GEO -- with user-settable thrust and
% user-settable initial/final orbits. Prints a Table-3-style row and
% (optionally) writes a static trajectory PNG and a burn/coast movie
% (mp4+gif). Modeled on PSR/run_psr.m and elfo/elfo_run_one.m: a single
% PARAMETERS block the user edits and runs (`run_gergaud`, no args), that
% is ALSO callable as `run_gergaud(opts)` (field-by-field override of the
% PARAMETERS block) so it is unit-testable without editing the file.
%
% This is a THIN FRONT DOOR onto the already-built MEE/L-domain thrust-
% ladder campaign (mee_seed.m / casadi_lt_mee.m / homotopy_mee.m /
% run_mintime_mee.m / run_transfer_mee.m / psr_mee_refine.m) -- it adds no
% new solver physics, just endpoint resolution, cache-vs-solve selection,
% and Table-3 row/plot/movie assembly.
%
% PARAMETERS block (edit directly for interactive use; every field below
% is also an `opts.<field>` override for `run_gergaud(opts)`):
%   thrustN     - max thrust [N]: one of 10 | 5 | 2.5 | 1 | 0.5 | 0.2 | 0.1
%                 (only 10/5/2.5/1/0.5 are campaign-certified; see the
%                 per-rung recipe map below)                    [default 10]
%   P0_km       - initial orbit semi-latus rectum P0 [km]   [default 11625]
%   e0          - initial eccentricity                       [default 0.75]
%   i0_deg      - initial inclination [deg]                      [default 7]
%   Pf_km       - final orbit P [km]                        [default 42165]
%   ef          - final eccentricity                              [default 0]
%   if_deg      - final inclination [deg]                          [default 0]
%   (P0_km,e0,i0_deg)=(11625,0.75,7) and (Pf_km,ef,if_deg)=(42165,0,0) are
%   the paper's own endpoints (GTO-like ellipse -> GEO); leaving both at
%   their defaults is what lets 'auto' mode reuse the certified caches.
%   ctf         - t_f / t_f,min ratio for the fixed-tf fuel solve [default 1.5]
%   nodesPerRev - collocation node density [nodes/rev]             [default 25]
%   maxIter     - IPOPT iteration cap per continuation step      [default 1500]
%   runMode     - 'auto' | 'solve' | 'probe' (see RUN MODES below) [default 'auto']
%   makeMovie   - write results/gergaud_<tag>.{mp4,gif}            [default true]
%   makePlot    - write results/gergaud_<tag>.png                  [default true]
%   m0kg        - initial spacecraft mass [kg]                    [default 1500]
%   ispS        - specific impulse [s]                            [default 2000]
%   returnOnly  - test hook: return the row struct and skip both
%                 plot and movie regardless of makePlot/makeMovie [default false]
%
% RUN MODES:
%   'auto'  (default) -- if BOTH endpoints are the paper defaults AND a
%           certified cache exists for thrustN, LOAD it and build the row
%           with NO solve. If either condition fails (custom endpoints, or
%           no cache for this thrustN), behaves exactly like 'solve'.
%   'solve' -- always runs the live pipeline (run_mintime_mee anchor ->
%           run_transfer_mee fixed-tf fuel homotopy -> psr_mee_refine for
%           thrustN<=1 N), ignoring any cache. Required for any custom
%           endpoint (no cache applies to a non-default target/initial orbit).
%   'probe' -- research mode: forces a live solve (like 'solve') and prints
%           an up-front WARNING that thrustN<0.5 N (0.2/0.1 N) was never
%           certified in this campaign before attempting it. Reports
%           row.certified honestly either way -- never fabricates a row.
%
% ENDPOINT RESOLUTION (default-preserving; see mee_seed.m/casadi_lt_mee.m):
%   initElems = [] when (P0_km,e0,i0_deg) match the paper defaults exactly
%   (byte-preserving legacy literal inside mee_seed.m); otherwise
%   initElems = [P0_km/LU; e0; 0; tan(deg2rad(i0_deg)/2); 0; 1; 0].
%   xf = [1;0;0;0;0] when (Pf_km,ef,if_deg) match GEO exactly; otherwise
%   xf = [Pf_km/LU; ef; 0; tan(deg2rad(if_deg)/2); 0]. LU = 42165 km always
%   (kepler_lt_params.m's fixed length unit; a custom final orbit is
%   expressed in that SAME unit, never rescaled).
%   isDefaultEndpoints = isempty(initElems) && isequal(xf,[1;0;0;0;0]).
%   SCOPE CAVEAT: the solver/seed were validated for GEO-like (near-
%   circular, near-equatorial) targets. A significantly eccentric/inclined/
%   retrograde custom final orbit is research-probe territory -- this
%   script reports whether the live solve certified rather than presuming
%   any custom target converges. Custom endpoints get a hashed tag suffix
%   (endpoint_hash_suffix, below) so a custom run's cache files never
%   collide with the default-endpoint certified caches.
%   PSR + CUSTOM ENDPOINTS: psr_mee_refine.m always re-solves against the
%   DEFAULT GEO terminal [1;0;0;0;0] (it has no xf passthrough). So for a
%   custom-endpoint run at thrustN<=1 N (where the default-endpoint recipe
%   would call psr_mee_refine), the PSR step is SKIPPED instead -- calling
%   it would silently re-terminate the trajectory at GEO while the row
%   stays labeled with the user's custom target. The reported row is the
%   un-refined run_transfer_mee fuel solve, which already terminates
%   correctly at the custom xf; row.note carries an explicit
%   "PSR switch-refinement skipped for custom endpoints" caveat. Default
%   endpoints are unaffected -- PSR still runs exactly as before.
%
% PER-RUNG RECIPE MAP (honest, encoded here; see README.md/DESIGN_thrust_
% ladder.md for the full campaign record; DEFAULT ENDPOINTS ONLY -- custom
% endpoints at thrustN<=1 N skip the PSR step, see PSR + CUSTOM ENDPOINTS
% above):
%   T [N]     | recipe                                    | status
%   10/5/2.5  | run_mintime_mee + run_transfer_mee         | clean, cached
%   1         | + psr_mee_refine switch-refinement         | headline is PSR round 2
%   0.5       | anchor-free R0-law tfmin (446.27 ND, the   | anchor footnoted an
%             | anchor solve hits a documented conditioning|  ESTIMATE, PSR is
%             | wall) + psr_mee_refine                     |  budget-limited
%   0.2/0.1   | live probe only (anchor + fuel + PSR,      | NEVER certified;
%             | same recipe, honestly attempted)           |  reports certified=false
%             |                                             |  rather than a
%             |                                             |  fabricated row
%
% OUTPUTS:
%   row - the gergaud_row() struct (see gergaud_row.m for all fields);
%         ALSO printed via gergaud_row_str() unless the printed row would
%         be redundant with returnOnly's test-only usage (it still prints;
%         returnOnly only suppresses viz). When ~row.certified, the printed
%         block carries an "UNCERTIFIED -- <note>" banner (gergaud_row_str).
%
% SIDE EFFECTS (unless returnOnly): writes results/gergaud_<tag>.png
% (makePlot) and results/gergaud_<tag>.{mp4,gif} (makeMovie), tag =
% mee_fuel_tag(thrustN) [+ endpoint hash suffix if custom].
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, "Low Thrust Minimum-Fuel Orbital
%       Transfer: A Homotopic Approach," JGCD 27(6), 2004, Table 3.
%   [2] earth_elliptic_to_geo/README.md, DESIGN_thrust_ladder.md (campaign
%       record + the six binding footnotes the recipe map above summarizes).
%   [3] earth_elliptic_to_geo/run_mintime_mee.m, run_transfer_mee.m,
%       psr_mee_refine.m (the three live-pipeline stages).
%   [4] earth_elliptic_to_geo/mee_res_to_cart_res.m, transfer_movie.m,
%       gergaud_plot.m (visualization adapter + renderers).
%   [5] earth_elliptic_to_geo/gergaud_row.m, gergaud_row_str.m (row
%       assembly/formatting, pure functions this script wraps).
%   [6] .superpowers/sdd/task-7-brief.md (this script's spec).

if nargin < 1, opts = struct(); end

%% =======================================================================
%% PARAMETERS  (edit this section for interactive use)
%% =======================================================================
thrustN     = 10;        % max thrust [N]
P0_km       = 11625;     % initial orbit P [km]              (paper default)
e0          = 0.75;      % initial eccentricity               (paper default)
i0_deg      = 7;         % initial inclination [deg]          (paper default)
Pf_km       = 42165;     % final orbit P [km]                 (GEO default)
ef          = 0;         % final eccentricity                 (GEO default)
if_deg      = 0;         % final inclination [deg]            (GEO default)
c_tf        = 1.5;       % t_f / t_f,min ratio
nodesPerRev = 25;        % collocation node density [nodes/rev]
maxIter     = 1500;      % IPOPT iteration cap per continuation step
runMode     = 'auto';    % 'auto' | 'solve' | 'probe'
makeMovie   = true;      % write results/gergaud_<tag>.{mp4,gif}
makePlot    = true;      % write results/gergaud_<tag>.png
m0kg        = 1500;      % initial spacecraft mass [kg]
ispS        = 2000;      % specific impulse [s]
returnOnly  = false;     % test hook: return row, skip plot+movie

%% ---- opts overrides (field-by-field; only used when called with opts) --
d = @(f,v) optdef(opts, f, v);
thrustN     = d('thrustN',     thrustN);
P0_km       = d('P0_km',       P0_km);
e0          = d('e0',          e0);
i0_deg      = d('i0_deg',      i0_deg);
Pf_km       = d('Pf_km',       Pf_km);
ef          = d('ef',          ef);
if_deg      = d('if_deg',      if_deg);
c_tf        = d('ctf',         c_tf);
nodesPerRev = d('nodesPerRev', nodesPerRev);
maxIter     = d('maxIter',     maxIter);
runMode     = d('runMode',     runMode);
makeMovie   = d('makeMovie',   makeMovie);
makePlot    = d('makePlot',    makePlot);
m0kg        = d('m0kg',        m0kg);
ispS        = d('ispS',        ispS);
returnOnly  = d('returnOnly',  returnOnly);

assert(any(strcmpi(runMode, {'auto','solve','probe'})), ...
    'run_gergaud:badRunMode', 'runMode must be ''auto''|''solve''|''probe'', got ''%s''', runMode);

resDir = fullfile(module_root(), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

%% =======================================================================
%% 1-2. Resolve endpoints (default-preserving; see header ENDPOINT RESOLUTION)
%% =======================================================================
par = kepler_lt_params(thrustN, m0kg, ispS);   % LU_km is fixed (42165), independent of thrustN
LU  = par.LU_km;

tol = 1e-9;
isPaperInit = abs(P0_km-11625) < tol && abs(e0-0.75) < tol && abs(i0_deg-7) < tol;
if isPaperInit
    initElems = [];                            % legacy literal inside mee_seed.m, byte-preserving
else
    initElems = [P0_km/LU; e0; 0; tan(deg2rad(i0_deg)/2); 0; 1; 0];
end

isGEOFinal = abs(Pf_km-42165) < tol && abs(ef-0) < tol && abs(if_deg-0) < tol;
if isGEOFinal
    xf = [1;0;0;0;0];
else
    xf = [Pf_km/LU; ef; 0; tan(deg2rad(if_deg)/2); 0];
end

isDefaultEndpoints = isempty(initElems) && isequal(xf, [1;0;0;0;0]);

if isDefaultEndpoints
    epSuffix = '';
else
    epSuffix = endpoint_hash_suffix(P0_km, e0, i0_deg, Pf_km, ef, if_deg);
end
fuelTag = [mee_fuel_tag(thrustN) epSuffix];

%% =======================================================================
%% 3. Mode: cache-hit (auto+default) vs live solve (auto-fallback/solve/probe)
%% =======================================================================
usedCache = false;
nrm = [];  %#ok<NASGU>
tfmin_ND = NaN; tfmin_h = NaN; anchorNote = '';

if strcmpi(runMode, 'auto') && isDefaultEndpoints
    [cacheFile, layout] = gergaud_cache_file(resDir, thrustN);
    if isfile(cacheFile)
        fprintf('[run_gergaud] AUTO: T=%g N -- loading certified cache %s (NO solve)\n', ...
            thrustN, cacheFile);
        nrm = load_cache_normalized(cacheFile, layout);
        [tfmin_ND, tfmin_h, anchorNote] = gergaud_cached_anchor(resDir, thrustN, m0kg, ispS);
        usedCache = true;
    else
        fprintf(['[run_gergaud] AUTO: no certified cache for T=%g N at the default endpoints ' ...
                 '-- falling back to a live SOLVE\n'], thrustN);
    end
end

if ~usedCache
    if strcmpi(runMode, 'probe')
        warning('run_gergaud:probeWall', ['PROBE MODE: forcing a live solve (T=%g N). 0.2/0.1 N ' ...
            'were NEVER certified in this campaign -- the 0.5 N min-time anchor already hit a ' ...
            'documented conditioning wall (README.md footnote 1) -- so any rung at or below that ' ...
            'is research-probe territory, not a reproduction of a known-good result. ' ...
            'row.certified will be reported honestly regardless of what this solve reaches.'], thrustN);
    end
    [nrm, tfmin_ND, tfmin_h, anchorNote] = solve_live(thrustN, c_tf, nodesPerRev, maxIter, ...
        m0kg, ispS, xf, initElems, fuelTag, isDefaultEndpoints); %#ok<ASGLU>
end

%% =======================================================================
%% 4-5. Assemble + print the Table-3 row
%% =======================================================================
noteParts = {};
if ~isempty(anchorNote),        noteParts{end+1} = anchorNote; end
if ~isempty(nrm.note),          noteParts{end+1} = nrm.note;   end
if ~isDefaultEndpoints,         noteParts{end+1} = 'custom endpoints (research-probe scope)'; end
note = strjoin(noteParts, '; ');

inp = struct('thrustN', thrustN, 'tfmin_ND', tfmin_ND, 'ctf', nrm.ctf, 'tf_ND', nrm.tf, ...
    'm_f_kg', nrm.m_f_kg, 'switches', nrm.switches, 'revs', nrm.revs, 'edge', nrm.edge, ...
    'incl_deg', nrm.incDeg, 'defect', nrm.defect, 'certified', logical(nrm.certified), ...
    'note', note, 'm0kg', m0kg, 'ispS', ispS);
row = gergaud_row(inp);
fprintf('%s', gergaud_row_str(row));

if returnOnly
    return;
end

%% =======================================================================
%% 6. Plot + movie (via the Cartesian adapter)
%% =======================================================================
outStem = fullfile(resDir, ['gergaud_' fuelTag]);
if makePlot || makeMovie
    if isempty(nrm.Xmee)
        fprintf(['[run_gergaud] no trajectory state available (uncertified/failed solve) -- ' ...
                 'skipping plot/movie\n']);
    else
        cartRes = mee_res_to_cart_res(nrm.Xmee, nrm.Umee, nrm.dL, nrm.sigma, thrustN, nrm.ctf, par.mu);
        if makePlot
            gergaud_plot(cartRes, [outStem '.png']);
        end
        if makeMovie
            transfer_movie(cartRes, outStem);
        end
    end
end

end

% =============================================================================
function suf = endpoint_hash_suffix(P0_km, e0, i0_deg, Pf_km, ef, if_deg)
% ENDPOINT_HASH_SUFFIX  Short deterministic tag suffix for a custom endpoint
% set, so a custom run's cache files (fuel/mintime/PSR) never collide with
% the default-endpoint certified caches under the same base tag.
vals = [P0_km, e0*1000, i0_deg*1000, Pf_km, ef*1000, if_deg*1000];
h = mod(round(abs(sum(vals .* (1:numel(vals)) * 97.531))), 1e6);
suf = sprintf('_ep%06d', h);
end

% =============================================================================
function [cacheFile, layout] = gergaud_cache_file(resDir, thrustN)
% GERGAUD_CACHE_FILE  Where the certified headline fuel result for thrustN
% lives, ONLY valid for the default (paper/GEO) endpoints. 1 N/0.5 N use the
% PSR-final layout (psr_mee_refine.m output); everything else uses the plain
% run_transfer_mee.m layout. thrustN values with no campaign cache (0.2/0.1
% N) simply resolve to a non-existent file -- the caller falls through to a
% live solve.
tag = mee_fuel_tag(thrustN);
if abs(thrustN - 1) < 1e-6 || abs(thrustN - 0.5) < 1e-6
    cacheFile = fullfile(resDir, [tag '_PSR_psr_final.mat']);
    layout    = 'psr';
else
    cacheFile = fullfile(resDir, [tag '.mat']);
    layout    = 'plain';
end
end

% =============================================================================
function nrm = load_cache_normalized(cacheFile, layout)
% LOAD_CACHE_NORMALIZED  Normalize either cache layout run_transfer_mee.m
% (plain: res.cfg/.fuel/.report/.sigma/.tf) or psr_mee_refine.m (psr:
% out.finalOut/.finalSigma/.certified) produces into one common struct:
% {Xmee,Umee,dL,sigma,m_f_kg,switches,revs,edge,incDeg,defect,certified,
% tf,ctf,note}. Verified against the actual cache files on disk (2026-07-18):
% MEE_M2_10N.mat/MEE_M2_5N.mat/MEE_M2_2p5N.mat (plain) and
% MEE_M2_1N_PSR_psr_final.mat/MEE_M2_0p5N_PSR_psr_final.mat (psr).
S = load(cacheFile);
switch layout
    case 'plain'
        res = S.res;
        nrm = struct('Xmee', res.fuel.X, 'Umee', res.fuel.U, 'dL', res.fuel.dL, ...
            'sigma', res.sigma(:), 'm_f_kg', res.report.m_f_kg, 'switches', res.report.switches, ...
            'revs', res.report.revs, 'edge', res.report.edge, 'incDeg', res.report.incDeg, ...
            'defect', res.report.defect, 'certified', logical(res.report.certified), ...
            'tf', res.tf, 'ctf', optdef(res.cfg, 'ctf', 1.5), 'note', '');
    case 'psr'
        out = S.out;
        fo  = out.finalOut;
        nrm = struct('Xmee', fo.X, 'Umee', fo.U, 'dL', fo.dL, 'sigma', out.finalSigma(:), ...
            'm_f_kg', fo.m_f_kg, 'switches', fo.switches, 'revs', fo.dL/(2*pi), 'edge', fo.edge, ...
            'incDeg', fo.incDeg, 'defect', fo.maxDefect, 'certified', logical(out.certified), ...
            'tf', fo.tf, 'ctf', 1.5, 'note', '');
    otherwise
        error('run_gergaud:badLayout', 'load_cache_normalized: unknown layout ''%s''', layout);
end
end

% =============================================================================
function [tfmin_ND, tfmin_h, note] = gergaud_cached_anchor(resDir, thrustN, m0kg, ispS)
% GERGAUD_CACHED_ANCHOR  Min-time anchor lookup for the printed row, default
% endpoints only. 0.5 N has NO certified anchor (documented conditioning
% wall, README.md footnote 1) -- returns the campaign's R0-law estimate
% (tfmin ~= 446.27 ND) with an honest note instead. Everything else reads
% results/MEE_mintime_T<round(10*thrustN)>.mat (out.tfmin/.tfmin_h).
note = '';
if abs(thrustN - 0.5) < 1e-6
    par = kepler_lt_params(thrustN, m0kg, ispS);
    tfmin_ND = 446.27;
    tfmin_h  = tfmin_ND * par.TU_s / 3600;
    note = '0.5 N: R0-law tfmin estimate (anchor-free)';
    return;
end
mtFile = fullfile(resDir, sprintf('MEE_mintime_T%d.mat', round(10*thrustN)));
if isfile(mtFile)
    S = load(mtFile);
    tfmin_ND = S.out.tfmin;
    tfmin_h  = S.out.tfmin_h;
else
    tfmin_ND = NaN; tfmin_h = NaN;
    note = sprintf('T=%g N: no cached min-time anchor found (%s)', thrustN, mtFile);
end
end

% =============================================================================
function [nrm, tfmin_ND, tfmin_h, note] = solve_live(thrustN, ctf, nodesPerRev, maxIter, ...
    m0kg, ispS, xf, initElems, fuelTag, isDefaultEndpoints)
% SOLVE_LIVE  The live pipeline: min-time anchor -> fixed-tf fuel homotopy
% -> (for thrustN<=1 N) PSR switch-aware mesh refinement. Honest on
% failure: any stage that throws or fails to certify yields a
% certified=false nrm carrying a descriptive note rather than propagating
% the exception or fabricating a row (gergaud_row_str already banners
% uncertified rows on its own).
%
% Per-rung recipe (see run_gergaud.m header PER-RUNG RECIPE MAP):
%   thrustN == 0.5  -> the live min-time anchor is a documented conditioning
%                      wall (README.md footnote 1); SKIP it and use the
%                      R0-law estimate tfmin ~= 446.27 ND directly.
%   otherwise       -> run_mintime_mee (works cleanly for 10/5/2.5/1 N;
%                      attempted honestly, may fail, for 0.2/0.1 N probes).
%   thrustN <= 1 N  -> psr_mee_refine on top of the fixed-tf fuel solve
%                      (matches the 1 N/0.5 N/deeper campaign recipe).
note = '';
mtTagBase = sprintf('MEE_mintime_T%d', round(10*thrustN));
if isDefaultEndpoints
    mtTag = mtTagBase;
else
    fuelBase = mee_fuel_tag(thrustN);
    epSuf    = fuelTag(numel(fuelBase)+1:end);   % '' or '_ep######'
    mtTag    = [mtTagBase epSuf];
end

isHalfN = abs(thrustN - 0.5) < 1e-6;
if thrustN < 0.5 - 1e-9
    fprintf(['[run_gergaud] NOTE: T=%g N was NEVER certified in this campaign (the 0.5 N ' ...
             'min-time anchor already hit a conditioning wall; README.md footnote 1) -- this ' ...
             'is a live research-probe solve; certified will be reported honestly.\n'], thrustN);
end

if isHalfN
    par0 = kepler_lt_params(thrustN, m0kg, ispS);
    tfmin_ND = 446.27;
    tfmin_h  = tfmin_ND * par0.TU_s / 3600;
    if isDefaultEndpoints
        note = '0.5 N: R0-law tfmin estimate (anchor-free)';
    else
        note = ['0.5 N: R0-law tfmin ESTIMATE carried over from the default-endpoint campaign ' ...
                 'fit -- NOT re-derived for these custom endpoints, treat as approximate'];
    end
    fprintf(['[run_gergaud] T=0.5 N: SKIPPING the live min-time anchor solve (documented ' ...
             'conditioning wall) -- using R0-law estimate tfmin=%.2f ND\n'], tfmin_ND);
else
    fprintf('[run_gergaud] SOLVE: T=%g N min-time anchor...\n', thrustN);
    try
        mtCfg = struct('m0kg', m0kg, 'ispS', ispS, 'maxIter', maxIter, 'tag', mtTag, ...
            'xf', xf, 'initElems', initElems);
        anchorOut = run_mintime_mee(thrustN, nodesPerRev, mtCfg);
        tfmin_ND = anchorOut.tfmin; tfmin_h = anchorOut.tfmin_h;
        if ~anchorOut.certified
            note = sprintf('T=%g N: min-time anchor did NOT certify', thrustN);
        end
    catch ME_anchor
        tfmin_ND = NaN; tfmin_h = NaN;
        note = sprintf(['T=%g N: min-time anchor solve THREW (%s) -- likely the known deep-' ...
            'ladder conditioning wall'], thrustN, ME_anchor.message);
        % note carried in the function's own note output only (see caller's
        % noteParts assembly) -- do NOT also fold it into nrm.note, or the
        % printed UNCERTIFIED banner would show the same diagnostic twice.
        nrm = uncertified_nrm(ctf, '');
        return;
    end
end

fprintf('[run_gergaud] SOLVE: T=%g N fuel homotopy (ctf=%.2f)...\n', thrustN, ctf);
try
    fuelCfg = struct('thrustN', thrustN, 'ctf', ctf, 'tfMinAnchor', tfmin_ND, 'tag', fuelTag, ...
        'nodesPerRev', nodesPerRev, 'maxIter', maxIter, 'm0kg', m0kg, 'ispS', ispS, ...
        'xf', xf, 'initElems', initElems);
    res = run_transfer_mee(fuelCfg);
catch ME_fuel
    note = strtrim([note ' ' sprintf('T=%g N: fuel homotopy THREW (%s)', thrustN, ME_fuel.message)]);
    % note carried in the function's own note output only -- see the
    % ME_anchor catch above for why nrm.note stays empty here.
    nrm = uncertified_nrm(ctf, '');
    return;
end

needsPSR = (thrustN <= 1 + 1e-9);
if needsPSR && ~isDefaultEndpoints
    % FIX I-1: psr_mee_refine.m (and its internal solve_psr_round) do NOT
    % accept/forward a terminal target xf -- every PSR round re-solves
    % against the hardcoded default GEO terminal [1;0;0;0;0]. Calling it
    % here for a custom-endpoint run would silently refine the trajectory
    % back toward GEO while the row/plot/movie stay labeled with the
    % user's custom target: a silent wrong answer. Rather than threading
    % xf into the validated PSR core, SKIP the PSR step for custom
    % endpoints and report the un-refined fuel solve, which already
    % terminates correctly at the custom xf (run_transfer_mee threads xf
    % through the fuel homotopy) -- just not switch-sharpened.
    fprintf(['[run_gergaud] T=%g N: PSR switch-refinement SKIPPED for custom endpoints -- ' ...
             'psr_mee_refine always re-targets the default GEO terminal, so refining here ' ...
             'would silently re-terminate at GEO instead of the requested custom xf. ' ...
             'Reporting the un-refined fuel solve (already correctly targeted).\n'], thrustN);
    nrm = res_to_nrm(res, ctf, ['PSR switch-refinement skipped for custom endpoints ' ...
        '(research-probe): reported solution is the un-refined fuel solve']);
elseif needsPSR
    fprintf('[run_gergaud] SOLVE: T=%g N -- applying PSR switch-refinement recipe...\n', thrustN);
    try
        psrOpts = struct('tag', [fuelTag '_PSR'], 'maxRounds', 4, 'maxIter', maxIter, 'nbr', 2);
        psrOut  = psr_mee_refine(res, psrOpts);
        fo = psrOut.finalOut;
        nrm = struct('Xmee', fo.X, 'Umee', fo.U, 'dL', fo.dL, 'sigma', psrOut.finalSigma(:), ...
            'm_f_kg', fo.m_f_kg, 'switches', fo.switches, 'revs', fo.dL/(2*pi), 'edge', fo.edge, ...
            'incDeg', fo.incDeg, 'defect', fo.maxDefect, 'certified', logical(psrOut.certified), ...
            'tf', fo.tf, 'ctf', ctf, 'note', '');
        if ~psrOut.certified
            nrm.note = sprintf('T=%g N: PSR stopped (%s) without certifying', thrustN, psrOut.stopReason);
        end
    catch ME_psr
        note = strtrim([note ' ' sprintf(['T=%g N: PSR refinement THREW (%s) -- falling back ' ...
            'to the pre-PSR fuel solve'], thrustN, ME_psr.message)]);
        % note carried in the function's own note output only -- see the
        % ME_anchor catch above for why nrm.note stays empty here.
        nrm = res_to_nrm(res, ctf, '');
    end
else
    % note (if any, e.g. an anchor-did-not-certify message from above) is
    % carried in the function's own note output only -- same reasoning as
    % the catch blocks above, so it isn't printed twice by the caller.
    nrm = res_to_nrm(res, ctf, '');
end
end

% =============================================================================
function nrm = res_to_nrm(res, ctf, note)
% RES_TO_NRM  Normalize a run_transfer_mee.m 'res' struct (plain layout)
% into the common {Xmee,Umee,dL,sigma,m_f_kg,switches,revs,edge,incDeg,
% defect,certified,tf,ctf,note} struct -- same field set load_cache_
% normalized produces, so downstream row/plot/movie code is layout-agnostic.
nrm = struct('Xmee', res.fuel.X, 'Umee', res.fuel.U, 'dL', res.fuel.dL, 'sigma', res.sigma(:), ...
    'm_f_kg', res.report.m_f_kg, 'switches', res.report.switches, 'revs', res.report.revs, ...
    'edge', res.report.edge, 'incDeg', res.report.incDeg, 'defect', res.report.defect, ...
    'certified', logical(res.report.certified), 'tf', res.tf, 'ctf', ctf, 'note', note);
end

% =============================================================================
function nrm = uncertified_nrm(ctf, note)
% UNCERTIFIED_NRM  Placeholder normalized result for a stage that threw
% before producing any trajectory -- certified=false, no state to plot/
% animate, but still a well-formed struct so gergaud_row/gergaud_row_str
% can print an honest UNCERTIFIED row instead of the caller crashing.
nrm = struct('Xmee', [], 'Umee', [], 'dL', NaN, 'sigma', [], 'm_f_kg', NaN, 'switches', NaN, ...
    'revs', NaN, 'edge', NaN, 'incDeg', NaN, 'defect', NaN, 'certified', false, 'tf', NaN, ...
    'ctf', ctf, 'note', note);
end
