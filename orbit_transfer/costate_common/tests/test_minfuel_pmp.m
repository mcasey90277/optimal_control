function ok = test_minfuel_pmp()
%% Purpose:
%
%   Tests cr3bp_minfuel_pmp -- the smoothed energy->fuel PMP field with two
%   smoothing families (Bertrand-Epenoy 'eps', PLQ Huber 'huber'):
%
%     1. eps = 1 reproduces cr3bp_minenergy_pmp (F, A, s) to round-off --
%        L_eps = (1-eps)s + eps s^2 equals s^2 there.
%     2. The throttle law is the true argmin of H(s) on [0,1] for BOTH
%        families (checked against a dense scan of H(s) at random states).
%     3. Both families approach the SAME bang-bang limit: at small
%        parameter, fields agree with each other away from the switch
%        (Q > 1.05 and Q < 0.95), where Q = T(|lam_v|/m + lam_m/c).
%     4. The AD Jacobian matches central finite differences of F.
%     5. The Huber law's PREDICTED discontinuity: s* jumps from ~kappa to 1
%        as Q crosses 1 (the affine-tail degeneracy; recorded theory, must
%        hold so the race measures the real family, not a bug).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All cases passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/01/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
mu_ = 0.012150585609624;  T_ = 0.1756418;  c_ = 8.673746;
rng(7);

% base catalog-like state; lam_v magnitude rescaled per-case to set Q:
yb = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 0.97; ...
      15.6; 32.9; -0.09; -0.10; 0.045; -0.00015; 0.13];
setQ = @(y, Q) setQfun(y, Q, T_, c_);

%% 1. eps = 1 == min-energy field:
yt = setQ(yb, 1.3);                       % interior throttle for energy
[Fe, Ae, axe] = cr3bp_minenergy_pmp(yt, T_, c_, mu_);
[Ff, Af, axf] = cr3bp_minfuel_pmp(yt, T_, c_, mu_, struct('family','eps','p',1));
ok = chk(ok, max(abs(Ff - Fe)) < 1e-12*(1 + max(abs(Fe))) && ...
             max(abs(Af(:) - Ae(:))) < 1e-11*(1 + max(abs(Ae(:)))) && ...
             abs(axf.s - axe.s) < 1e-13, ...
         sprintf('eps=1 reproduces min-energy (dF %.1e, ds %.1e)', ...
                 max(abs(Ff - Fe)), abs(axf.s - axe.s)));

%% 2. throttle law = argmin_s H(s) (dense scan), both families:
fams = { struct('family','eps',  'p',0.35); ...
         struct('family','eps',  'p',0.05); ...
         struct('family','huber','p',0.6);  ...
         struct('family','huber','p',0.15) };
worst = 0;
for kf = 1:numel(fams)
    for Qt = [0.3, 0.8, 0.999, 1.2, 2.5]
        y = setQ(yb, Qt);
        [~, ~, ax] = cr3bp_minfuel_pmp(y, T_, c_, mu_, fams{kf});
        sg = linspace(0, 1, 20001);
        Hg = Lof(sg, fams{kf}) - sg*ax.Q;         % H(s) modulo s-free terms
        [~, ib] = min(Hg);
        % distance in H-VALUE (flat segments make s itself ambiguous):
        Hstar = Lof(ax.s, fams{kf}) - ax.s*ax.Q;
        worst = max(worst, Hstar - Hg(ib));
    end
end
ok = chk(ok, worst < 1e-7, ...
         sprintf('throttle law is the H(s) argmin, both families (worst dH %.1e)', worst));

%% 3. common bang-bang limit away from the switch. Above the switch the
%  agreement is exact; BELOW it huber converges only FIRST-ORDER in kappa
%  (its law s* = kappa*Q never coasts exactly -- no dead zone, a real
%  family property the race must know), so the tolerance there is O(p):
p3 = 1e-6;
for Qt = [0.7, 1.4]
    y = setQ(yb, Qt);
    [F1] = cr3bp_minfuel_pmp(y, T_, c_, mu_, struct('family','eps','p',p3));
    [F2] = cr3bp_minfuel_pmp(y, T_, c_, mu_, struct('family','huber','p',p3));
    tol = 10*p3*(1 + max(abs(F1)));       % first-order-in-kappa envelope
    ok = chk(ok, max(abs(F1 - F2)) < tol, ...
             sprintf('eps/huber share the bang limit at Q = %.1f (dF %.1e < %.0e)', ...
                     Qt, max(abs(F1 - F2)), tol));
end

%% 4. AD Jacobian vs central finite differences:
for kf = [1 3]
    y = setQ(yb, 1.25);
    [~, A] = cr3bp_minfuel_pmp(y, T_, c_, mu_, fams{kf});
    Afd = zeros(14);
    h = 1e-6;
    for kc = 1:14
        yp = y; ym = y; yp(kc) = yp(kc) + h; ym(kc) = ym(kc) - h;
        Afd(:,kc) = (cr3bp_minfuel_pmp(yp, T_, c_, mu_, fams{kf}) - ...
                     cr3bp_minfuel_pmp(ym, T_, c_, mu_, fams{kf}))/(2*h);
    end
    dA = max(abs(A(:) - Afd(:)))/(1 + max(abs(A(:))));
    ok = chk(ok, dA < 1e-5, ...
             sprintf('%s AD Jacobian vs FD (rel %.1e)', fams{kf}.family, dA));
end

%% 5. the predicted Huber jump at Q = 1 (theory check, kappa = 0.3):
kap = 0.3;
[~,~,axLo] = cr3bp_minfuel_pmp(setQ(yb, 0.995), T_, c_, mu_, ...
                               struct('family','huber','p',kap));
[~,~,axHi] = cr3bp_minfuel_pmp(setQ(yb, 1.005), T_, c_, mu_, ...
                               struct('family','huber','p',kap));
ok = chk(ok, abs(axLo.s - kap*0.995) < 1e-3 && abs(axHi.s - 1) < 1e-12, ...
         sprintf('huber law jumps kappa->1 across Q=1 (s %.3f -> %.3f)', ...
                 axLo.s, axHi.s));

if ok, fprintf('TEST_MINFUEL_PMP: ALL PASS\n');
else,  fprintf('TEST_MINFUEL_PMP: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label.
% OUTPUTS: ok updated.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end

% ------------------------------------------------------------------------
function y = setQfun(y, Q, T, c)
% SETQFUN  Rescale lam_v so the switch quantity Q = T(|lam_v|/m + lam_m/c)
% hits the requested value.  INPUTS: y [14x1]; Q; T; c.  OUTPUTS: y.
m = y(7);  lm = y(14);
nlvWant = m*(Q/T - lm/c);
assert(nlvWant > 0, 'setQfun: requested Q unreachable with this lam_m');
y(11:13) = y(11:13) * nlvWant/sqrt(sum(y(11:13).^2));
end

% ------------------------------------------------------------------------
function L = Lof(s, fam)
% LOF  The running cost L(s) of a smoothing family.
% INPUTS: s [1xN]; fam struct('family','p').  OUTPUTS: L [1xN].
p = fam.p;
switch fam.family
    case 'eps'
        L = (1 - p)*s + p*s.^2;
    case 'huber'
        L = (s <= p).*(s.^2/(2*p)) + (s > p).*(s - p/2);
end
end
