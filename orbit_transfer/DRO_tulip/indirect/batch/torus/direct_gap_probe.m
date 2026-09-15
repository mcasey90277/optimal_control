%% DIRECT_GAP_PROBE  Branch-blind direct solves (Hermite-Simpson + Sundman,
% IPOPT) at the five uncertified arrival phases of the 24x24 library, each
% warm-started from the CERTIFIED indirect solution on either side of the gap
% (col 18 at sA 0.7837, 24.74 d; col 21 at 0.9087, 26.43 d), with the 1900 km
% lunar clearance enforced as a path constraint. Two control targets (the
% seeds' own phases) validate the seeding. Progress: results/direct_gap_probe.txt
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
D = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/direct';
addpath(ind, D, fullfile(D, 'lib'), fullfile(D, 'certify'), '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
logF = fullfile(D, 'results', 'direct_gap_probe.txt');
lg = @(varargin) logline(logF, sprintf(varargin{:}));
S = load(fullfile(ind, 'results_fine', 'arrival_sheet_70mN_nA24.mat'));  S = S.S;
[B, ~] = arclength_arrival('setup', struct('sD', 0));
mu = B.mu;  Tmax = B.Tnd;  c = B.cnd;  lStar = 389703.264829278;  tStar = 382981.289129055;  rMoonKm = 1737.4;
rv0 = B.stateD(0);  rv0 = rv0(1:6);
N = 800;  floorKm = 1900 - rMoonKm;                 % the certifier's clearance, as altitude
seedCols = [18 21];
targets = [0.7837 0.9087 0.8254 0.8671 0.9504 0.9921 0.0337];
seeds = struct('col', {}, 'sA', {}, 'X', {}, 'U', {}, 'tf', {});
for j = seedCols
    cnd_ = S.cand{j};  k = find([cnd_.ok], 1);  z = cnd_(k).z;
    [tauR, rvR] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(:); 1; z(1:7)], Tmax, c, mu);
    sN = linspace(0, tauR(end), N+1);
    sd.col = j;  sd.sA = S.sA(j);
    sd.X = interp1(tauR, rvR(:, 1:7), sN, 'spline').';
    LV = interp1(tauR, rvR(:, 11:13), sN, 'spline');
    sd.U = [(-LV./max(vecnorm(LV, 2, 2), eps)).'; ones(1, N+1)];
    sd.tf = z(8);
    seeds(end+1) = sd; %#ok<SAGROW>
    lg('seed from col %d (sA %.4f): t_f %.3f d', j, S.sA(j), z(8)*tStar/86400);
end
lg('DIRECT GAP PROBE: N = %d, minAlt %.1f km, %d targets x %d seeds', N, floorKm, numel(targets), numel(seeds));
P = struct('sA', {}, 'seedCol', {}, 'success', {}, 'tfDays', {}, 'perisKm', {}, 'globKm', {}, 'globMs', {}, 'maxDefect', {}, 'wall', {}, 'o', {});
for sA = targets
    rvf = B.stateA(mod(sA, 1));  rvf = rvf(1:6);
    for s = seeds
        t0 = tic;  rec = struct('sA', sA, 'seedCol', s.col, 'success', false, 'tfDays', NaN, 'perisKm', NaN, ...
                                'globKm', Inf, 'globMs', Inf, 'maxDefect', NaN, 'wall', NaN, 'o', []);
        try
            o = casadi_mintime_dro(rv0, rvf, Tmax, c, mu, N, s.X, s.U, s.tf, struct('maxIter', 3000, ...
                    'scheme', 'hermite-simpson', 'sundman', true, 'returnModel', true, 'minAltKm', floorKm, 'maxCpuSec', 900));
            rec.success = o.success;  rec.maxDefect = o.maxDefect;  rec.tfDays = o.tf*tStar/86400;
            rec.perisKm = min(vecnorm(o.X(1:3, :) - [1-mu; 0; 0], 2, 1))*lStar - rMoonKm;
            if o.success && o.maxDefect < 1e-9
                C = certify_dro_mintime(o, struct('muStar', mu, 'lStar', lStar, 'tStar', tStar), Tmax, c, ...
                                        struct('tfRef', [], 'verbose', false, 'posTolKm', inf));
                rec.globKm = C.globKm;  rec.globMs = C.gates(strcmp({C.gates.id}, 'G1bv')).value;
            end
            o = rmfield(o, intersect(fieldnames(o), {'model'}));  rec.o = o;
        catch ME
            lg('  sA %.4f seed col %d: ERROR %s', sA, s.col, ME.message);
        end
        rec.wall = toc(t0);
        P(end+1) = rec; %#ok<SAGROW>
        lg('  sA %.4f  seed col %2d:  success=%d  t_f %8.3f d  peris %7.0f km  glob %.3f km / %.3f m/s  defect %.1e  (%.0f s)', ...
           sA, s.col, rec.success, rec.tfDays, rec.perisKm, rec.globKm, rec.globMs, rec.maxDefect, rec.wall);
        save(fullfile(D, 'results', 'direct_gap_probe.mat'), 'P', 'targets', 'seedCols', 'N', 'floorKm');
    end
end
lg('DIRECT GAP PROBE DONE');
function logline(f, s)
fid = fopen(f, 'a');  fprintf(fid, '%s %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), s);  fclose(fid);
fprintf('%s\n', s);
end
