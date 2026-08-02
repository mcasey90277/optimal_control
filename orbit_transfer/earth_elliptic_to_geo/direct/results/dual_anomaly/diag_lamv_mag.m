% DIAG_LAMV_MAG  Are the lamV duals tiny (direction ill-conditioned) on burns?
% Compare ||lamV|| against the stationarity coupling scale and against angle.
cd('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo');
S = load('results/M1_3d_fixedLf.mat');  o = S.res.fuel;
N = size(o.lamDef,2);  ss = o.U(4,1:N);  burn = find(ss > 0.5);
lamV = o.lamDef(4:6,:);  lVn = sqrt(sum(lamV.^2,1));
lamAllMag = sqrt(sum(o.lamDef.^2,1));            % full 9-row dual magnitude
cosb = zeros(1,N);
for k = 1:N, cosb(k) = dot(o.U(1:3,k), -lamV(:,k))/max(lVn(k),1e-30); end
if mean(cosb(burn)) < 0, cosb = -cosb; end
ang = real(acosd(max(-1,min(1,cosb(burn)))));
fprintf('||lamV|| on burns: min=%.2e med=%.2e max=%.2e\n', ...
        min(lVn(burn)), median(lVn(burn)), max(lVn(burn)));
fprintf('||lamAll|| on burns: med=%.2e   ratio lamV/lamAll med=%.2e\n', ...
        median(lamAllMag(burn)), median(lVn(burn)./lamAllMag(burn)));
% correlation: are the worst angles at the smallest ||lamV||?
[~, isrt] = sort(lVn(burn));
angSorted = ang(isrt);
q = round(numel(isrt)/4);
fprintf('mean angle by ||lamV|| quartile (small->large): %.1f %.1f %.1f %.1f deg\n', ...
        mean(angSorted(1:q)), mean(angSorted(q+1:2*q)), ...
        mean(angSorted(2*q+1:3*q)), mean(angSorted(3*q+1:end)));
