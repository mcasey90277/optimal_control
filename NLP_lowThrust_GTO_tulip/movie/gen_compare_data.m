% gen_compare_data  Assemble the three solutions (min-time, min-energy,
% min-fuel) into one .mat for the energy-solo and three-way movies.

sp  = [fileparts(mfilename('fullpath')) filesep];
mf  = load([sp 'minfuel_movie_data.mat']);       % min-time (yF) + min-fuel (out) + tulip
en  = load([sp 'energy_pipeline.mat']);          % N=4000 direct min-energy solve

muStar = mf.muStar;  tStar = mf.tStar;

% --- min-time: full spiral, throttle always 1 ------------------------------
rMT = mf.yF(:,1:3).';                 % 3 x M
tMT = mf.tauF(:).';
sMT = ones(1, numel(tMT));
tfMT = mf.tfMinFull;

% --- min-energy: full spiral, smooth throttle ramp -------------------------
rME = en.out.X(1:3,:);
sME = en.out.U(4,:);
tME = en.out.tauMesh(:).';
tfME = en.out.tf_ND;

% --- min-fuel: arrival leg only (bang-bang) --------------------------------
rMF = mf.out.X(1:3,:);
sMF = mf.out.U(4,:);
tMF = mf.out.tauMesh(:).';
tfMF = tMF(end);
legStart = mf.legStart;               % where the min-fuel leg begins (tau)

% --- fixed scene geometry --------------------------------------------------
earth = [-muStar, 0, 0];
moon  = [1 - muStar, 0, 0];
yTul  = mf.yTul;                       % tulip reference curve

save([sp 'compare_data.mat'], 'rMT','tMT','sMT','tfMT', ...
     'rME','sME','tME','tfME', 'rMF','sMF','tMF','tfMF','legStart', ...
     'earth','moon','yTul','muStar','tStar');
fprintf('SAVED compare_data.mat\n');
fprintf('min-time  tf=%.3f d,  min-energy tf=%.3f d,  min-fuel leg tf=%.3f d (from tau=%.1f)\n', ...
        tfMT*tStar/86400, tfME*tStar/86400, tfMF*tStar/86400, legStart);
