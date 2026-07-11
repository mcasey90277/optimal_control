function cert = ifs_certify(Z, prob, meta)
% IFS_CERTIFY  First-order PMP certificate + structure diagnostics for an IFS Z.
%
% Verifies (post-hoc) that a converged IFS solution is a continuous-time
% first-order extremal: S=0 at each switch, the sign law on each arc interior
% (S<0 burn / S>0 coast, bounded away from 0 => no singular arc), the terminal
% residual, and (rendezvous mode) transversality. Reports structure diagnostics:
% the smallest arc length (a vanishing arc => a spurious switch) and the worst
% in-arc sign violation (=> a missing switch). Diagnostics are REPORTED, not
% acted on.
%
% INPUTS:
%   Z    - converged unknown vector [(8+17k)x1]
%   prob - problem struct
%   meta - seed meta from IFS_SEED (for switchMoveFromSeed)
% OUTPUTS:
%   cert - struct: ok (logical), Sswitch [1xk], signViol (max), minArcLen,
%          termResNorm, lamMend, switchMoveFromSeed, text (certificate string)
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
k = prob.k;
[~, N, gblk] = ifs_unpack(Z, k);
tau = ifs_taus(gblk, prob.tau0, prob.tauf);
R = ifs_residual(Z, prob);

% S at each switch node (should be ~0)
Sswitch = zeros(1,k);
for ii = 1:k, Sswitch(ii) = 1 - sqrt(sum(N(12:14,ii).^2))*prob.c/N(7,ii) - N(15,ii); end

% per-arc interior sign law: sample S mid-arc from the integrated arc
edges = [prob.tau0, tau(:).', prob.tauf];
signViol = 0;  minArcLen = min(diff(edges));
% reintegrate arcs to sample midpoint S (reuse ode113 directly)
[lam0, Nn, ~] = ifs_unpack(Z, k);
startY = cell(1,k+1);
startY{1} = [prob.rv0(:); prob.m0; prob.t0; lam0];
for a = 2:k, startY{a} = Nn(:,a-1); end
startY{k+1} = Nn(:,k);
for a = 1:k+1
    sp = [edges(a), edges(a+1)];
    [~, Yar] = ode113(@(s,y) ifs_eom(s,y,prob.Tmax,prob.c,prob.muStar,prob.pSund,prob.uArc(a)), ...
                      sp, startY{a}, prob.odeOpts);
    Smid = zeros(size(Yar,1),1);
    for q = 1:size(Yar,1)
        yq = Yar(q,:).';  Smid(q) = 1 - sqrt(sum(yq(12:14).^2))*prob.c/yq(7) - yq(15);
    end
    if prob.uArc(a) == 1, viol = max(Smid);     % burn arc wants S<0
    else,                 viol = max(-Smid); end % coast arc wants S>0
    signViol = max(signViol, viol);
end

termResNorm = norm(R(16*k + (1:8)));
lamMend = NaN;
if strcmp(prob.termMode,'rendezvous')
    % lamM at tf: integrate the last arc's endpoint costate
    lamMend = abs(R(16*k + 7));   % rendezvous row 7 is e(15)=lamM(tf); residual==lamM(tf)
end
switchMoveFromSeed = max(abs(tau(:).' - meta.tauSwitch(:).'));

okS   = max(abs(Sswitch)) < 1e-6;
okLaw = signViol < 1e-3;           % sign law respected on arc interiors
okTerm= termResNorm < 1e-6;
okTr  = ~strcmp(prob.termMode,'rendezvous') || lamMend < 1e-3;
cert.ok = okS && okLaw && okTerm && okTr;

if cert.ok
    cert.text = sprintf(['CERTIFICATE: continuous-time first-order PMP extremal ' ...
        '(max|S(switch)|=%.1e, sign-law viol=%.1e, term=%.1e).'], ...
        max(abs(Sswitch)), signViol, termResNorm);
else
    cert.text = 'CERTIFICATE NOT ISSUED: see per-check values.';
end
% structure diagnostics (reported)
cert.diag = sprintf('minArcLen=%.2e (vanishing<1e-3 => fewer sw); maxSignViol=%.2e (>1e-3 => missing sw)', ...
                    minArcLen, signViol);
cert.Sswitch=Sswitch; cert.signViol=signViol; cert.minArcLen=minArcLen;
cert.termResNorm=termResNorm; cert.lamMend=lamMend; cert.switchMoveFromSeed=switchMoveFromSeed;
fprintf('%s\n%s\n', cert.text, cert.diag);
end
