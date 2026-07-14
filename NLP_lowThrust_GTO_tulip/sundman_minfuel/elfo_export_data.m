function dataFile = elfo_export_data(solFile, dataDir, opts)
% ELFO_EXPORT_DATA  Data products for a GTO->ELFO (free-t_f two-primary) min-fuel
% solution -- the ELFO analog of PSR/psr_export_data.
%
% Writes ONE self-contained data-product .mat that mirrors psr_export_data's
% layout (mesh / traj / ctrl / costate / pmp / scal / const / provenance), with
% two model differences that MATTER:
%   (1) Costates come from the SOLVER'S OWN KKT duals (out.lamDef, the two-primary
%       free-t_f discrete costates), NOT the tulip single-primary fixed-t_f
%       dual->costate map (sms_seed_duals) -- that map assumes a different model
%       and would give wrong costates here.
%   (2) The PMP certificate is the solver's scale-INVARIANT primer alignment
%       (mean angle between the NLP thrust and the costate primer -lamV/||lamV||
%       on burn arcs; ~0 = optimal) plus the terminal mass-costate transversality
%       proxy. The switching function S is a BEST-EFFORT single-parameter (beta)
%       scaling of the raw costates against the throttle -- honestly labeled, NOT
%       the tulip's adjudicated dual-S certified count.
%
% INPUTS:
%   solFile - a minfuel_<target>_*.mat from gen_elfo_minfuel: needs top-level
%             out (with X[9xnN],U[4xnN],lamDef[9xN],primerAlignDeg,lamMassEnd),
%             sigma, tauf0, rv0, rvf, tf, moonZone, pSund, qSund, target, factor.
%   dataDir - destination dir (created if missing) [default ../PSR_data]
%   opts    - (optional): .quiet [false]
%
% OUTPUTS:
%   dataFile - written path psr_data_<target>_tf<f>_sw<k>_minEps<e>.mat.
%
% REFERENCES:
%   [1] PSR/psr_export_data.m (the tulip analog this mirrors).
%   [2] casadi_energy_freetf.m (the two-primary solver whose costates are used).

if nargin < 2 || isempty(dataDir)
    dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'PSR_data');
end
if nargin < 3, opts = struct(); end
quiet = isfield(opts,'quiet') && opts.quiet;
if ~exist(dataDir,'dir'), mkdir(dataDir); end

here = fileparts(mfilename('fullpath'));
p = cr3bp_lt_params(0.025, 15, 2100);
lStar = p.lStar;  tStar = p.tStar;  m0kg = p.m0kg;  cEx = p.c;  Tmax = p.Tmax;

S = load(solFile);
out = S.out;  X = out.X;  U = out.U;  nN = size(X,2);
sigma = S.sigma(:).';  tauf0 = S.tauf0;  rv0 = S.rv0;  rvf = S.rvf;  rvfC = rvf(:);
target = S.target;  factor = S.factor;  epsMin = S.epsilon;  tf = S.tf;
moonZone = S.moonZone;  pSund = S.pSund;  qSund = S.qSund;

r = X(1:3,:);  v = X(4:6,:);  m = X(7,:);  t = X(8,:);  cScale = X(9,:);
al = U(1:3,:);  s = U(4,:);
tau = sigma*tauf0;  tDays = t*tStar/86400;

% ---- mesh / traj / ctrl -----------------------------------------------------
mesh = struct('sigma',sigma,'tau',tau,'tauf0',tauf0,'t',t,'tDays',tDays, ...
              'pSund',pSund,'qSund',qSund,'moonZone',moonZone,'nN',nN);
traj = struct('r',r,'v',v,'m',m,'cScale',cScale,'X',X);

% throttle switch times (model-independent): tau where s crosses 0.5
sb = s > 0.5;  cr = find(diff(sb) ~= 0);
tauSwitch = zeros(1,numel(cr));
for k = 1:numel(cr)
    j = cr(k);  a0 = s(j)-0.5;  a1 = s(j+1)-0.5;  w = a0/(a0-a1);
    tauSwitch(k) = tau(j) + w*(tau(j+1)-tau(j));
end
ctrl = struct('alpha',al,'s',s,'nSwitchThrottle',numel(tauSwitch),'tauSwitchThrottle',tauSwitch);

% ---- costates from the solver's own two-primary KKT duals -------------------
% map interval duals lamDef [9xN] -> node costates [9xnN] (adjacent-interval avg)
costate = struct();  pmpBeta = NaN;  signAgree = NaN;
if isfield(out,'lamDef') && ~isempty(out.lamDef) && size(out.lamDef,1) == 9
    lamD = out.lamDef;  Nint = size(lamD,2);
    lamNode = zeros(9, nN);
    lamNode(:,1) = lamD(:,1);  lamNode(:,end) = lamD(:,end);
    lamNode(:,2:Nint) = 0.5*(lamD(:,1:Nint-1) + lamD(:,2:Nint));
    lamV = lamNode(4:6,:);  lamM = lamNode(7,:);
    normLamV = sqrt(sum(lamV.^2,1));
    % best-effort switching function S = 1 - beta*q, q = (Tmax/m)||lamV|| + (Tmax/cEx)lamM.
    % Fit beta (and a global sign) to agree with the throttle (S<0 burn, S>0 coast).
    q = (Tmax./m).*normLamV + (Tmax/cEx).*lamM;
    want = double(s <= 0.5) - double(s > 0.5);        % +1 coast, -1 burn = target sign(S)
    bestAgree = -inf;  bestBeta = NaN;
    for sg = [1 -1]
        qs = sg*q;  qpos = qs(qs>0);
        if isempty(qpos), continue; end
        for bb = 1./quantile(qpos, linspace(0.05,0.95,25))
            agree = mean(sign(1 - bb*qs) == want);
            if agree > bestAgree, bestAgree = agree; bestBeta = sg*bb; end
        end
    end
    pmpBeta = bestBeta;  signAgree = bestAgree;
    S_switch = 1 - bestBeta*q;
    costate = struct('lam', lamNode, 'lamInterval', lamD, 'S', S_switch, ...
        'beta', pmpBeta, 'signAgreePct', 100*signAgree, ...
        'primerAlignDeg', out.primerAlignDeg, 'lamMassEnd', out.lamMassEnd, ...
        'note', ['solver two-primary KKT costates (up to +mesh-scale/global-sign); ' ...
                 'S is a best-effort 1-param beta fit, NOT the tulip adjudicated dual-S']);
elseif ~quiet
    warning('elfo_export_data:noCostates','solution out.lamDef absent/wrong size -- costate layer skipped');
end

% ---- first-order PMP diagnostics -------------------------------------------
mf = m(end);
dVtot = cEx*log(1/mf)*lStar/tStar;              % rocket-equation dV (km/s)
pmp = struct('lamMassEnd', getfielddef(out,'lamMassEnd',NaN), ...
             'primerAlignDeg', getfielddef(out,'primerAlignDeg',NaN), ...
             'termPosErr', norm(X(1:3,end)-rvfC(1:3)), ...
             'termVelErr', norm(X(4:6,end)-rvfC(4:6)), ...
             'termTimeErr', t(end)-S.tf, 'SsignAgreePct', 100*signAgree);

% ---- scalars / constants / provenance --------------------------------------
scal = struct('target',target,'factor',factor,'tf',S.tf,'tf_days',S.tf*tStar/86400, ...
              'epsMin',epsMin,'dV',dVtot,'prop_kg',m0kg*(1-mf),'mf',mf, ...
              'maxDefect',getfielddef(out,'maxDefect',NaN),'switches',numel(tauSwitch), ...
              'edge',mean(s>0.95|s<0.05));
const = p;  const.tfMin = 6.2906939607;  const.moonZone = moonZone;  const.pSund = pSund;  const.qSund = qSund;
provenance = struct('date',datestr(now,'yyyy-mm-dd HH:MM:SS'), ... %#ok<TNOW1,DATST>
    'source',char(solFile),'gitHash',git_hash_local(here), ...
    'solver','casadi_energy_freetf (free-t_f, two-primary Sundman)', ...
    'pipeline','run_elfo_minfuel');

% ---- write ------------------------------------------------------------------
fTag = strrep(sprintf('%.3f', factor), '.', 'p');
eTag = strrep(sprintf('%g', epsMin), '.', 'p');
dataFile = fullfile(dataDir, sprintf('psr_data_%s_tf%s_sw%d_minEps%s.mat', ...
                    target, fTag, numel(tauSwitch), eTag));
save(dataFile, 'out','sigma','tauf0','rv0','rvf','factor','target','tf', ...  %#ok<USENS> seed layer
     'mesh','traj','ctrl','costate','pmp','scal','const','provenance');
if ~quiet
    fprintf('WROTE %s\n', dataFile);
    fprintf('  %s tf=%.2f d  dV=%.4f km/s  prop=%.4f kg  sw=%d  edge=%.1f%%\n', ...
        target, scal.tf_days, scal.dV, scal.prop_kg, scal.switches, 100*scal.edge);
    fprintf('  PMP: primerAlign=%.3f deg  |lamM(end)|=%.2g  S-sign agree=%.1f%% (beta=%.4g)\n', ...
        pmp.primerAlignDeg, abs(pmp.lamMassEnd), 100*signAgree, pmpBeta);
end
end

% ===========================================================================
function v = getfielddef(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function h = git_hash_local(dir0)
[rc, so] = system(sprintf('cd "%s" && git rev-parse --short HEAD 2>/dev/null', dir0));
if rc == 0, h = strtrim(so); else, h = 'unknown'; end
end
