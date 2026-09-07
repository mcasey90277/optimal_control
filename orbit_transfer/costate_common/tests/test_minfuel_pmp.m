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
         struct('family','huber','p',0.15); ...
         struct('family','huberc','p',0.6);                 % delta defaults to p
         struct('family','huberc','p',0.15,'delta',0.05) };
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
         sprintf('throttle law is the H(s) argmin, all families (worst dH %.1e)', worst));

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
for kf = [1 3 5]
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

%% 6. the Huber-eps HYBRID 'huberc' (FINDINGS 24 cure): Huber's core
%  s = p*Q below Q = 1, then a CONTINUOUS ramp p -> 1 over Q in [1, 1+delta]
%  instead of the jump. (a) continuous at both knees; (b) delta -> 0
%  recovers huber away from the ramp; (c) the ramp really is the law.
hc = struct('family','huberc','p',0.3,'delta',0.1);
sAt = @(Qt, fam) sAtFun(setQ(yb, Qt), T_, c_, mu_, fam);
ok = chk(ok, abs(sAt(1-1e-7, hc) - sAt(1+1e-7, hc)) < 1e-5 && ...
             abs(sAt(1.1-1e-7, hc) - sAt(1.1+1e-7, hc)) < 1e-5, ...
         sprintf('huberc continuous at Q=1 and Q=1+delta (jumps %.1e, %.1e)', ...
                 abs(sAt(1-1e-7,hc)-sAt(1+1e-7,hc)), abs(sAt(1.1-1e-7,hc)-sAt(1.1+1e-7,hc))));
hc0 = struct('family','huberc','p',0.3,'delta',1e-9);
hu  = struct('family','huber','p',0.3);
ok = chk(ok, abs(sAt(0.8, hc0) - sAt(0.8, hu)) < 1e-12 && abs(sAt(1.5, hc0) - sAt(1.5, hu)) < 1e-12, ...
         'huberc(delta->0) == huber away from the ramp');
ok = chk(ok, abs(sAt(1.05, hc) - (0.3 + 0.7*0.5)) < 1e-9, ...
         sprintf('huberc ramp: s(Q=1+delta/2) = p + (1-p)/2 = %.4f', sAt(1.05, hc)));

%% 7. Astra review #2 (2026-09-06): huberc validation, p = 1 dispatch, sRaw
% (a) delta must be a finite positive scalar:
bad = 0;
for d = {0, -0.1, NaN, Inf}
    try, cr3bp_minfuel_pmp(yb, T_, c_, mu_, struct('family','huberc','p',0.3,'delta',d{1})); catch, bad = bad + 1; end
end
ok = chk(ok, bad == 4, sprintf('huberc rejects delta in {0, <0, NaN, Inf} (%d/4 rejected)', bad));
% (b) p = 1: huberc must equal huber(p = 1) exactly (no 1/(1-p) evaluation):
y = setQ(yb, 0.6);
[F1, A1] = cr3bp_minfuel_pmp(y, T_, c_, mu_, struct('family','huberc','p',1));
[F2, A2] = cr3bp_minfuel_pmp(y, T_, c_, mu_, struct('family','huber','p',1));
ok = chk(ok, all(isfinite(F1)) && max(abs(F1-F2)) < 1e-13 && max(abs(A1(:)-A2(:))) < 1e-12, ...
         sprintf('huberc(p=1) == huber(p=1): dF %.1e dA %.1e', max(abs(F1-F2)), max(abs(A1(:)-A2(:)))));
% (c) sRaw on the ramp is the piecewise inverse law (unclipped), not p*Q:
[~, ~, axr] = cr3bp_minfuel_pmp(setQ(yb, 1.15), T_, c_, mu_, struct('family','huberc','p',0.3,'delta',0.3));
ok = chk(ok, abs(axr.sRaw - 0.65) < 1e-9, sprintf('huberc sRaw on the ramp = p + (1-p)(Q-1)/delta = %.4f', axr.sRaw));

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
    case 'huberc'
        d = p;  if isfield(fam, 'delta') && ~isempty(fam.delta), d = fam.delta; end
        L = (s <= p).*(s.^2/(2*p)) + ...
            (s > p).*(p/2 + (s - p) + d*(s - p).^2/(2*(1 - p)));
end
end

function s = sAtFun(y, T, c, mu, fam)
% SATFUN  Throttle of a smoothing family at a state.
% INPUTS: y [14x1]; T; c; mu; fam struct.  OUTPUTS: s double.
[~, ~, ax] = cr3bp_minfuel_pmp(y, T, c, mu, fam);
s = ax.s;
end
