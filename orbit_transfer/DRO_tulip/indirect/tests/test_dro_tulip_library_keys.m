function ok = test_dro_tulip_library_keys()
% TEST_DRO_TULIP_LIBRARY_KEYS  every entry of dro_tulip_library() arrives
% where its (sD, sA) key says it does. Two entries were once loaded from
% files solved before the 2026-09-12 phase-to-state rule (periodic spline):
% their costates solved the files' own endpoints and missed the keyed
% arrival by 49 km and 22 km under the current rule. Each entry is flown from
% B.stateD(sD) and its miss at B.stateA(sA) must be < 1 km and < 0.1 m/s.
% Also: with opts.includeCatalog no two entries share a phase pair.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));                     % IND
cc = fullfile(fileparts(fileparts(here)), 'costate_common');
addpath(here, cc);
[B, ~] = arclength_arrival('setup');
phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'muStar', B.mu, 'lStar', 389703.264829278, ...
              'tStar', 382981.289129055, 'm0kg', 150);

lib = dro_tulip_library();
for k = 1:numel(lib)
    E = lib(k);
    rv0 = B.stateD(E.sD);  rvf = B.stateA(E.sA);
    F = fly_transfer(E.z(:), rv0(1:6), rvf(1:6), phys);
    rMiss = F.flyKm;   vMiss = F.flyVms;                                % km, m/s
    ok = chk(ok, rMiss < 1 && vMiss < 0.1, ...
        sprintf('%-22s (%.4f, %.4f): miss %.3f km / %.4f m/s', E.src, E.sD, E.sA, rMiss, vMiss));
end

lc = dro_tulip_library(here, struct('includeCatalog', true));
key = round([lc.sD; lc.sA]' * 1e9);
ok = chk(ok, size(unique(key, 'rows'), 1) == numel(lc), ...
    sprintf('includeCatalog: %d entries, all phase pairs distinct', numel(lc)));
nBad = 0;  worst = 0;
for k = 1:numel(lc)
    E = lc(k);
    rv0 = B.stateD(E.sD);  rvf = B.stateA(E.sA);
    F = fly_transfer(E.z(:), rv0(1:6), rvf(1:6), phys);
    worst = max(worst, F.flyKm);
    if ~(F.flyKm < 1 && F.flyVms < 0.1)
        nBad = nBad + 1;
        fprintf('    off key: %s (%.4f, %.4f): %.3f km / %.4f m/s\n', E.src, E.sD, E.sA, F.flyKm, F.flyVms);
    end
end
ok = chk(ok, nBad == 0, sprintf('includeCatalog: all %d entries arrive on their keys (worst %.3f km)', numel(lc), worst));
fprintf('test_dro_tulip_library_keys: %s\n', ternary(ok));
end

function ok = chk(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end

function s = ternary(ok)
if ok, s = 'ALL PASS'; else, s = 'FAILED'; end
end
