function results = thrust_continuation_minfuel_indirect(saveFile)
% THRUST_CONTINUATION_MINFUEL_INDIRECT  Zhang-style thrust homotopy for the
% min-fuel GTO -> tulip transfer in the Earth-Moon CR3BP.
%
% The fresh attack (after Zhang, Topputo, Bernelli-Zazzera & Zhao, JGCD 2015):
% start at HIGH thrust, where the transfer is only a few revolutions and the
% shooting basin is large, then step T_max DOWN, warm-starting the COSTATES
% from the previous thrust level. Costate continuation naturally grows the
% revolution/switch count as thrust drops -- the thing a direct trajectory
% warm start cannot do. This is the ingredient the earlier campaign never
% tried (we did tf- and energy->fuel continuation at FIXED thrust).
%
% Pipeline:
%   Phase 1  min-time up-continuation (25 mN -> high) for tf_min(T) + a
%            costate backbone at the top thrust (min-time is robust).
%   Phase 2  min-fuel at the top thrust, energy->fuel homotopy from a
%            rescaled min-time seed (few revs -> converges).
%   Phase 3  min-fuel down-continuation to 25 mN, warm-starting costates.
%
% Reuses SOLVE_TFMIN_INDIRECT, SOLVE_MINFUEL_INDIRECT, LT_PMP_EOM(_MINFUEL).
%
% INPUTS:
%   saveFile - (optional) .mat for incremental results
%
% OUTPUTS:
%   results - struct array per thrust level: .factor .Tmax_mN .tf .tf_min
%             .resNorm .switches .mProp_kg .dV_kms .lam0 .converged

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(saveFile)
    saveFile = fullfile(here, 'thrust_continuation_results.mat');
end
run(fullfile(here, 'setup_paths.m'));      % pumpkyn

muStar = 0.012150585609624;  lStar = 389703.264829278;  tStar = 382981.289129055;
m0kg = 15;  g0 = 9.80665*tStar^2/(1000*lStar);  c = (2100/tStar)*g0;
Tmax25 = (0.025/m0kg)*tStar^2/(lStar*1000);    % ND accel at nominal 25 mN

% --- endpoints (GTO -> tulip max-ydot point), same as run_gto_tulip_indirect
muEarth = 6.67384e-20*(1-muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;  ecc = (35786-350)/(2*sma);
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0 = pumpkyn.cr3bp.fromPCI(0, [r0, v0], muStar, tStar, lStar, 1);
[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
[~, rvTgt]   = pumpkyn.cr3bp.prop((5/6)*2*pi, x0Tulip, muStar);
[~, idxF]    = max(rvTgt(:,5));
rvf = rvTgt(idxF, :);

zGuess25 = [190.476497248065; -79.7064866984696; -0.430399154713168; ...
             0.301159446575878; 0.586671892449694; -0.00711582435720301; ...
             4.32931089137559; 6.29081541876621];

factors = fliplr(1.08.^(0:18));                % ~8% steps, 4.0x -> 1.0x (25 mN)
ctf     = 1.15;                                 % tf = ctf * tf_min at each level
optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% ================= Phase 1: min-time up-continuation ======================
upF = sort(factors);                            % ascending
tfMin = zeros(1, numel(factors));  costMT = zeros(7, numel(factors));
zmt = zGuess25;
fprintf('=== Phase 1: min-time thrust continuation (up) ===\n');
for iu = 1:numel(upF)
    f = upF(iu);  Tmax = f*Tmax25;
    [zmt, rn] = solve_tfmin_indirect(rv0, rvf, zmt, Tmax, c, muStar);
    j = find(abs(factors - f) < 1e-9, 1);
    tfMin(j) = zmt(8);  costMT(:,j) = zmt(1:7);
    fprintf('  f=%4.1f (%.0f mN): tf_min=%.4f ND, ||R||=%.2g\n', ...
            f, 25*f, zmt(8), rn);
end

% ================= Phases 2-3: min-fuel down-continuation ==================
results = struct('factor',{},'Tmax_mN',{},'tf',{},'tf_min',{},'resNorm',{}, ...
                 'switches',{},'mProp_kg',{},'dV_kms',{},'lam0',{},'converged',{});
lamPrev = [];
fprintf('\n=== Phases 2-3: min-fuel thrust continuation (down) ===\n');
for i = 1:numel(factors)
    f = factors(i);  Tmax = f*Tmax25;  tf = ctf*tfMin(i);
    if isempty(lamPrev)
        % bootstrap: rescale the min-time costates so S sits in the sensitive
        % zone (||lam_v|| c ~ 1 near switches), then full eps-homotopy
        lamMT = costMT(:,i);
        lamSeed = lamMT / (sqrt(sum(lamMT(4:6).^2))*c);
        eps0 = [1 0.3 0.1 0.03 0.01 3e-3 1e-3];
    else
        lamSeed = lamPrev;                       % costate continuation
        eps0 = [0.03 0.01 3e-3 1e-3];
    end
    ok = false;
    try
        [lamSol, rn] = solve_minfuel_indirect(rv0, 1, rvf, tf, lamSeed, ...
                            Tmax, c, muStar, eps0);
        % integrate, count switches, accounting
        [~, rvI] = ode113(@lt_pmp_eom_minfuel, [0 tf], [rv0(:); 1; lamSol], ...
                          optsInt, Tmax, c, muStar, 1e-3);
        lamvMag = sqrt(sum(rvI(:,11:13).^2, 2));
        S  = 1 - lamvMag.*c./rvI(:,7) - rvI(:,14);
        sw = sum(abs(diff(S > 0)));
        conv = rn < 1e-3;
        results(end+1) = struct('factor',f,'Tmax_mN',25*f,'tf',tf,'tf_min',tfMin(i), ...
            'resNorm',rn,'switches',sw,'mProp_kg',m0kg*(1-rvI(end,7)), ...
            'dV_kms',c*log(1/rvI(end,7))*lStar/tStar,'lam0',lamSol,'converged',conv); %#ok<AGROW>
        fprintf('  f=%4.1f (%.0f mN): tf=%.3f  ||R||=%.2g  switches=%3d  prop=%.4f kg  %s\n', ...
                f, 25*f, tf, rn, sw, m0kg*(1-rvI(end,7)), ternary(conv,'OK','(loose)'));
        if conv, lamPrev = lamSol;  ok = true; end   % advance only on convergence
    catch meErr
        fprintf('  f=%4.1f FAILED: %s\n', f, meErr.message);
    end
    save(saveFile, 'results', 'factors', 'tfMin');
    if ~ok && ~isempty(lamPrev)
        fprintf('    (kept previous costate seed)\n');
    end
end

fprintf('\n=== thrust-continuation summary (switches vs thrust) ===\n');
for k = 1:numel(results)
    fprintf('  %.0f mN: switches=%3d  ||R||=%.2g  prop=%.4f kg  %s\n', ...
            results(k).Tmax_mN, results(k).switches, results(k).resNorm, ...
            results(k).mProp_kg, ternary(results(k).converged,'','LOOSE'));
end
end

function v = ternary(cond, a, b)
if cond, v = a; else, v = b; end
end
