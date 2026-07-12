function [Z, prob, meta] = ifs_seed(matFile, opts)
% IFS_SEED  Build an IFS unknown vector + problem from a direct/PSR solution.
%
% Recovers node costates from the direct solution's KKT duals (mode-'d' map,
% SMS_SEED_DUALS), forms the switching function S(tau) on the node grid, takes
% its zero crossings as the switch times, samples the augmented state+costate at
% each switch as the node unknowns, and reads each arc's throttle from sign(S).
% 'full' poses the whole transfer (rendezvous terminal); 'window' extracts a
% one-switch sub-arc with fixed end states (fixedState terminal) for a
% ground-truth test.
%
% INPUTS:
%   matFile - direct solution .mat (out.X/out.U/out.lamDef, factor, tauf0, sigma)
%   opts    - struct: mode 'full'|'window' [required]; M dual-map arcs [40];
%             winSwitch switch index to center a window on [3]; winPad node pad
%             each side of the window switch [60]
% OUTPUTS:
%   Z    - seed unknown vector [(8+17k)x1]
%   prob - problem struct (see plan Shared data layout)
%   meta - struct: k, tauSwitch [1xk], uArc [1x(k+1)], seedResNorm, beta
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
if ~isfield(opts,'M'),        opts.M = 40;        end
if ~isfield(opts,'winSwitch'),opts.winSwitch = 3; end
if ~isfield(opts,'winPad'),   opts.winPad = 60;   end

[~, sd, info] = sms_seed_duals(matFile, opts.M, 1e-4, 'd');
Y16 = info.Y16;  tauN = info.tauN;  nN = size(Y16,2);
c = sd.c;

% switching function on the node grid and its zero crossings
Sn = 1 - sqrt(sum(Y16(12:14,:).^2,1)).*c./Y16(7,:) - Y16(15,:);
cr = find(diff(sign(Sn)) ~= 0);                 % node index before each crossing
tauCr = zeros(1, numel(cr));
for q = 1:numel(cr)
    kk = cr(q);
    tauCr(q) = tauN(kk) + (0-Sn(kk))*(tauN(kk+1)-tauN(kk))/(Sn(kk+1)-Sn(kk));
end

odeOpts = odeset('RelTol',1e-13,'AbsTol',1e-15);
sampleY = @(tt) interp1(tauN.', Y16.', tt(:), 'linear').';   % [16 x numel(tt)]
arcThrottle = @(tstart,tend) double( ...
    interp1(tauN, Sn, 0.5*(tstart+tend), 'linear') < 0 );     % burn where S<0

if strcmp(opts.mode, 'window')
    ws = opts.winSwitch;
    a0 = max(1, cr(ws) - opts.winPad);
    a1 = min(nN, cr(ws) + opts.winPad + 1);
    % ensure exactly ONE crossing inside [tauN(a0), tauN(a1)]
    inWin = tauCr > tauN(a0) & tauCr < tauN(a1);
    assert(nnz(inWin) == 1, 'window must contain exactly one switch (got %d); adjust winPad/winSwitch', nnz(inWin));
    Y0 = Y16(:, a0);  Yend = Y16(:, a1);
    tau1 = tauCr(inWin);
    N1 = sampleY(tau1);
    k = 1;
    prob = struct('rv0',Y0(1:6),'m0',Y0(7),'t0',Y0(8),'tau0',tauN(a0), ...
        'Tmax',sd.Tmax,'c',c,'muStar',sd.muStar,'pSund',sd.pSund, ...
        'tauf',tauN(a1),'k',1, ...
        'uArc',[arcThrottle(tauN(a0),tau1), arcThrottle(tau1,tauN(a1))], ...
        'termMode','fixedState','termTarget',Yend(1:8),'odeOpts',odeOpts);
    Z = ifs_pack(Y0(9:16), N1, ifs_gseed(tau1, tauN(a0), tauN(a1)));
    tauSwitch = tau1;
else   % full
    k = numel(tauCr);
    Y0 = Y16(:, 1);
    N  = sampleY(tauCr);                         % [16 x k]
    edges = [tauN(1), tauCr, tauN(end)];
    uArc = zeros(1, k+1);
    for a = 1:k+1, uArc(a) = arcThrottle(edges(a), edges(a+1)); end
    prob = struct('rv0',sd.rv0(:),'m0',1,'t0',0,'tau0',tauN(1), ...
        'Tmax',sd.Tmax,'c',c,'muStar',sd.muStar,'pSund',sd.pSund, ...
        'tauf',tauN(end),'k',k,'uArc',uArc, ...
        'termMode','rendezvous','rvf',sd.rvf(:),'tf',sd.tf,'odeOpts',odeOpts);
    Z = ifs_pack(Y0(9:16), N, ifs_gseed(tauCr(:), tauN(1), tauN(end)));
    tauSwitch = tauCr;
end

meta = struct('k',prob.k,'tauSwitch',tauSwitch,'uArc',prob.uArc,'beta',info.beta);
meta.seedResNorm = norm(ifs_residual(Z, prob));
end
