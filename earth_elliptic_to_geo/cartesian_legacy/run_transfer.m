function res = run_transfer(cfg)
% RUN_TRANSFER  One full pipeline: mintime anchor -> seed -> homotopy -> report.
%
% Stages: (1) cached free-L min-time at cfg.thrustN -> tfmin, dL_mt;
% (2) tf = ctf*tfmin; L_f = pi + (1.12*ctf+0.09)*dL_mt (paper law R2);
% (3) seed: tangential sbar=1/ctf bisected on L_f ('fixed' term) or plain
%     ('manifold'); or warm-start from cfg.seedMat with the light schedule;
% (4) homotopy eps 1->0; (5) structure report.
%
% INPUTS:  cfg - .thrustN .ctf .hx0 (0|0.0612) .term ('fixed'|'manifold')
%          .N (default 600) .tag (results filename stem) .seedMat (optional:
%          warm-start the homotopy from a prior result's res.fuel instead of
%          building a seed; uses the light schedule
%          [0.05 0.02 0.008 0.003 0.001 0]) .ispS (default 2000)
% OUTPUTS: res - .cfg .mintime .tf .Lf .fuel .tbl .report (ALWAYS returned;
%          saved to results/<tag>.mat ONLY when best.certified -- campaign
%          rule "never cache uncertified" (Task 14 controller triage): an
%          uncertified point warns and leaves no file, so a per-point resume
%          scan will retry it rather than silently adopting a bad optimum as
%          a downstream seed); .report = .revs .switches .m_f_kg .dV_kms
%          .edge .apoBurnRatio
%
% REFERENCES: [1] DESIGN.md secs 4-5.
resDir = fullfile(module_root(), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
d = @(f,v) optdef(cfg, f, v);
N = d('N', 600);  ispS = d('ispS', 2000);  seedMat = d('seedMat', '');

p  = kepler_lt_params(cfg.thrustN, 1500, ispS);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, cfg.hx0, 0, pi, p.mu);
rv0 = [r0; v0];
assert(ispS == 2000, 'run_transfer:ispMismatch', ...
    ['cfg.ispS=%g but run_mintime anchors are built/cached at Isp=2000 s only ' ...
     '(cache tag does not encode Isp) -- thread ispS through run_mintime before ' ...
     'using a different value'], ispS);
mt  = run_mintime(cfg.thrustN, cfg.hx0, N);
tf  = cfg.ctf * mt.tfmin;
Lf  = pi + (1.12*cfg.ctf + 0.09) * mt.dL_mt;
switch cfg.term
    case 'fixed',    term = geo_terminal('fixed', p, Lf);
    case 'manifold', term = geo_terminal('manifold', p, []);
end
ho = struct('par', p, 'rv0', rv0, 'maxIter', 1500);
if ~isempty(seedMat)                       % neighbor-style warm start
    if isfield(cfg,'N') && ~isempty(cfg.N)
        warning('run_transfer:seedMatMesh', ...
            'cfg.N ignored: mesh is inherited verbatim from seedMat (no-resample rule)');
    end
    S = load(seedMat);
    sg = S.res.sg;  tauf0 = S.res.fuel.tauf0;
    Xk = S.res.fuel.X;  Uk = S.res.fuel.U;
    Xk(8,:) = Xk(8,:) * (tf / Xk(8,end));  % rescale carried time if tf differs
    ho.sched = [0.05 0.02 0.008 0.003 0.001 0];
else
    so = struct('sbar', 1/cfg.ctf, 'tDur', tf, 'N', N);
    if strcmp(cfg.term, 'fixed'), so.targetLf = Lf; end
    [sg, Xk, Uk, tauf0] = seed_2body(p, rv0, so);
end
fprintf('RUN_TRANSFER %s: T=%g N, ctf=%.2f, tf=%.3f ND (%.1f h), Lf=%.2f rad\n', ...
        cfg.tag, cfg.thrustN, cfg.ctf, tf, tf*p.TU_s/3600, Lf);
[best, tbl] = homotopy_2body(sg, Xk, Uk, tauf0, term, tf, ho);

% structure report
Lun  = unwrap(atan2(best.X(2,:), best.X(1,:)));
revs = (Lun(end) - Lun(1)) / (2*pi);
rr   = sqrt(sum(best.X(1:3,:).^2, 1));
ss   = best.U(4,:);
nEarly = round(0.8 * numel(ss));           % exclude near-circular endgame
bMask  = ss(1:nEarly) > 0.5;
apoBurnRatio = median(rr(bMask)) / median(rr(~bMask));
report = struct('revs', revs, 'switches', best.switches, 'm_f_kg', best.m_f_kg, ...
    'dV_kms', best.dV_kms, 'edge', best.edge, 'apoBurnRatio', apoBurnRatio, ...
    'defect', best.maxDefect, 'certified', best.certified);
res = struct('cfg', cfg, 'mintime', mt, 'tf', tf, 'Lf', Lf, 'fuel', best, ...
             'tbl', tbl, 'report', report, 'sg', sg, 'rv0', rv0);
if best.certified
    save(fullfile(resDir, [cfg.tag '.mat']), 'res');
else
    warning('run_transfer:uncertified', ['%s: NOT saved (certified=0, ' ...
        'defect=%.2e) -- campaign rule: never cache uncertified results'], ...
        cfg.tag, best.maxDefect);
end
fprintf(['DONE %s: certified=%d revs=%.2f sw=%d edge=%.1f%% mf=%.2f kg ' ...
         'dV=%.3f km/s apoBurn=%.2f\n'], cfg.tag, report.certified, revs, ...
         best.switches, 100*best.edge, best.m_f_kg, best.dV_kms, apoBurnRatio);
end
