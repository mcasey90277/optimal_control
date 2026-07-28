% Discriminator: does the x8 18-switch structure survive on the COARSE mesh?
% If yes -> the certified row was in a worse BASIN, not under-resolved.
% If no  -> the 18-switch structure requires the fine mesh (a resolution effect).
D = '/Users/msc/Desktop/optimal_control/orbit_transfer/earth_elliptic_to_geo/direct';
cd(D); setup_paths();
addpath('/Users/msc/Desktop/optimal_control/orbit_transfer/verify_common');
setup_verify_common();
R = fullfile(D,'results','mesh_study');

F = load(fullfile(R,'MESH_MEE_M2_10N_x8.mat'));  fine = F.lvl;      % 18 sw
C = load(fullfile(R,'MESH_MEE_M2_10N_x1.mat'));  crse = C.lvl;      % 19 sw
saved = sosc_load_row(fullfile(D,'results','MEE_M2_10N.mat'));
par = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);

fprintf('fine  x8: N=%d  m_f=%.6f kg  sw=%d\n', fine.N, fine.out.m_f_kg, fine.out.switches);
fprintf('coarse x1: N=%d  m_f=%.6f kg  sw=%d\n', crse.N, crse.out.m_f_kg, crse.out.switches);

% project the FINE solution down onto the COARSE grid
sgC = crse.sigma(:);  sgF = fine.sigma(:);
Xd = interp1(sgF, fine.out.X.', sgC, 'pchip').';
Ud = interp1(sgF, fine.out.U.', sgC, 'pchip').';
n  = vecnorm(Ud(1:3,:),2,1);  Ud(1:3,:) = Ud(1:3,:)./max(n,eps);
Ud(4,:) = min(max(Ud(4,:),0),1);

sopts = struct('par',par,'mode','fixedtf','eps',0,'tfTarget',saved.tfTarget, ...
    'x0',saved.X(:,1),'xf',saved.xf,'maxIter',1500,'warmTight',true, ...
    'printLevel',0,'returnModel',false);
t0 = tic;
o = casadi_lt_mee(sgC, Xd, Ud, fine.out.dL, sopts);
fprintf('\nDOWN-PROJECTED to x1 (%.1f s): success=%d  %s\n', toc(t0), o.success, o.ipoptStatus);
fprintf('  m_f = %.6f kg   sw = %d   defect = %.2e   dV = %.6f\n', ...
        o.m_f_kg, o.switches, o.maxDefect, o.dV_kms);
fprintf('\n--- comparison on the SAME (production) mesh, N=%d ---\n', crse.N);
fprintf('  certified / warm-from-row seed : m_f = %.6f  sw = %d  dV = %.6f\n', ...
        crse.out.m_f_kg, crse.out.switches, crse.out.dV_kms);
fprintf('  seeded from the x8 solution    : m_f = %.6f  sw = %d  dV = %.6f\n', ...
        o.m_f_kg, o.switches, o.dV_kms);
d = o.m_f_kg - crse.out.m_f_kg;
fprintf('  difference                     : %+.6f kg  (%+.3e relative)\n', d, d/crse.out.m_f_kg);
if o.m_f_kg > crse.out.m_f_kg + 1e-4
    fprintf('\n=> BASIN: the better structure exists at PRODUCTION resolution.\n');
    fprintf('   The certified row is not under-resolved; it is in a worse basin.\n');
else
    fprintf('\n=> RESOLUTION: the coarse mesh cannot hold the fine structure.\n');
end
