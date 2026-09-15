%% DIRECT_GAP_PROBE4  The 4-day step at the fourth branch's lower end: the
% branch folds at 0.59, so cols 13/12/11 (0.5754/0.5337/0.4921) still carry
% the fast family at 21.5/20.8/20.2 d while col 14 (0.6171) has 17.27 d.
% Direct HS+Sundman solves at those phases, warm-started from the round-5
% sheet's certified roots at 0.6171 and 0.6587, with the clearance enforced.
% Solutions kept; certification is the sheet's job (mintime_70mN_direct_certified).
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
D = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/direct';
addpath(ind, D, fullfile(D, 'lib'), fullfile(D, 'certify'), '/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
logF = fullfile(D, 'results', 'direct_gap_probe4.txt');
lg = @(varargin) logline(logF, sprintf(varargin{:}));
S = load(fullfile(ind, 'results_fine_v5', 'arrival_sheet_70mN_nA24.mat'));  S = S.S;
[B, ~] = arclength_arrival('setup', struct('sD', 0));
mu = B.mu;  Tmax = B.Tnd;  c = B.cnd;  lStar = 389703.264829278;  tStar = 382981.289129055;  rMoonKm = 1737.4;
rv0 = B.stateD(0);  rv0 = rv0(1:6);
N = 800;  floorKm = 1900 - rMoonKm;
seeds = struct('sA', {}, 'X', {}, 'U', {}, 'tf', {});
for sAseed = [0.6171 0.6587]
    j = find(abs(mod(S.sA - sAseed + 0.5, 1) - 0.5) < 1e-4, 1);  cnd = S.cand{j};
    k = find([cnd.ok] & abs([cnd.tfDays] - S.TF(j)) < 1e-9, 1);  assert(~isempty(k), 'no certified winner at %.4f', sAseed);
    z = cnd(k).z;
    [tauR, rvR] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(:); 1; z(1:7)], Tmax, c, mu);
    sN = linspace(0, tauR(end), N+1);
    sd.sA = S.sA(j);  sd.X = interp1(tauR, rvR(:, 1:7), sN, 'spline').';
    LV = interp1(tauR, rvR(:, 11:13), sN, 'spline');
    sd.U = [(-LV./max(vecnorm(LV, 2, 2), eps)).'; ones(1, N+1)];  sd.tf = z(8);
    seeds(end+1) = sd;
    lg('seed: round-5 certified root at sA %.4f, t_f %.3f d', S.sA(j), z(8)*tStar/86400);
end
P = struct('sA', {}, 'seedSA', {}, 'success', {}, 'tfDays', {}, 'perisKm', {}, 'maxDefect', {}, 'wall', {}, 'o', {});
for sA = [0.5754 0.5337 0.4921]
    rvf = B.stateA(mod(sA, 1));  rvf = rvf(1:6);
    for s = seeds
        t0 = tic;  rec = struct('sA', sA, 'seedSA', s.sA, 'success', false, 'tfDays', NaN, 'perisKm', NaN, 'maxDefect', NaN, 'wall', NaN, 'o', []);
        try
            o = casadi_mintime_dro(rv0, rvf, Tmax, c, mu, N, s.X, s.U, s.tf, struct('maxIter', 3000, ...
                    'scheme', 'hermite-simpson', 'sundman', true, 'returnModel', true, 'minAltKm', floorKm, 'maxCpuSec', 900));
            if isfield(o, 'model'), o = rmfield(o, 'model'); end
            rec.o = o;  rec.success = o.success;  rec.maxDefect = o.maxDefect;  rec.tfDays = o.tf*tStar/86400;
            rec.perisKm = min(vecnorm(o.X(1:3, :) - [1-mu; 0; 0], 2, 1))*lStar - rMoonKm;
        catch ME
            lg('  sA %.4f seed %.4f: ERROR %s', sA, s.sA, ME.message);
        end
        rec.wall = toc(t0);  P(end+1) = rec;
        lg('  sA %.4f  seed %.4f:  success=%d  t_f %8.3f d  peris %7.0f km  defect %.1e  (%.0f s)', sA, s.sA, rec.success, rec.tfDays, rec.perisKm, rec.maxDefect, rec.wall);
        save(fullfile(D, 'results', 'direct_gap_probe4.mat'), 'P', 'N', 'floorKm');
    end
end
lg('DIRECT GAP PROBE4 DONE');
function logline(f, s)
fid = fopen(f, 'a');  fprintf(fid, '%s %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), s);  fclose(fid);  fprintf('%s\n', s);
end
