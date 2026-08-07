function B = survey_family_bounds(family, paramGrid, outMat)
%% Purpose:
%
%   Finds the "REASONABLE" members of ANY orbit family, by Darin's criteria:
%   periselene altitude >= 500 km (no lunar impact, with the same margin the
%   transfer campaigns use) and the whole orbit within 100 Mm of the Moon
%   (lunar vicinity, far enough out to take in L1/L2). Pure propagation via
%   get_family_orbit -- no optimization.
%
%   This generalizes the DRO/tulip survey that set the catalog's admissible
%   box; run it once per new family before building that family's sheets.
%
%% Inputs:
%
%  family                   char                    Family name understood
%                                                   by get_family_orbit
%
%  paramGrid                struct array            One element per candidate
%                                                   member: the parameter
%                                                   struct get_family_orbit
%                                                   expects (e.g. for halo:
%                                                   .tau .Lpt .pm)
%
%  outMat                   char                    Optional .mat for the
%                                                   survey table ('' = none)
%
%% Outputs:
%
%  B                        struct
%   .family                 char
%   .rows                   struct array             per member: params,
%                                                    periodDays, periAltKm,
%                                                    maxMoonKm, admissible
%   .admissible             logical vector
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, outMat = ''; end
   muStar = 0.012150585609624;
    lStar = 389703.264829278;
    tStar = 382981.289129055;
  rMoonKm = 1737.4;
altFloorKm = 500;
 maxDistKm = 100000;
    moonND = [1-muStar, 0, 0];

fprintf('%s family survey (%d candidates):\n', family, numel(paramGrid));
fprintf('   params                      period(d)   periAlt(km)   maxMoon(km)   admissible\n');
rows = [];
for k = 1:numel(paramGrid)
    p = paramGrid(k);
    try
        [tt, rv, info] = get_family_orbit(family, p);
        % PERIODICITY GUARD: an interpolated family seed that cont_np could
        % not truly converge produces a non-closing "orbit" whose metrics
        % are garbage (measured: periselene below the surface, distances of
        % 1e10 km). Label it, never trust it.
        perErr = norm(rv(end,1:6) - rv(1,1:6));
        % NOTE the negated form: a NaN closure (integrator blow-up mid-arc)
        % must fail this test too, and NaN passes any '>' comparison
        if ~(perErr < 1e-6) || any(~isfinite(rv(:)))
            fprintf('  %-26s  NON-PERIODIC (|closure| %.1e) -- excluded\n', ...
                    pstr(p), perErr);
            continue
        end
        dM = vecnorm(rv(:,1:3) - moonND, 2, 2) * lStar;
        pa = min(dM) - rMoonKm;
        mx = max(dM);
        ok = pa >= altFloorKm && mx <= maxDistKm;
        R = struct('params',p, 'periodDays',tt(end)*tStar/86400, ...
            'periAltKm',pa, 'maxMoonKm',mx, 'admissible',ok);
        if isempty(rows), rows = R; else, rows(end+1,1) = R; end %#ok<AGROW>
        fprintf('  %-26s  %8.3f   %11.0f   %11.0f   %s\n', pstr(p), ...
                tt(end)*tStar/86400, pa, mx, admstr(ok));
    catch e
        % never swallow the identity of an error -- a mislabeled catch cost
        % a debugging round (the real error was a struct-append mistake HERE,
        % reported as 'family generator failed')
        fprintf('  %-26s  FAILED: %s\n', pstr(p), e.message);
    end
end
adm = false(numel(rows),1);
for k = 1:numel(rows), adm(k) = rows(k).admissible; end
B = struct('family',family, 'rows',rows, 'admissible',adm, ...
           'criteria',struct('altFloorKm',altFloorKm,'maxDistKm',maxDistKm));
fprintf('ADMISSIBLE: %d of %d members\n', nnz(adm), numel(rows));
if ~isempty(outMat), save(outMat, '-struct', 'B'); end
end

% ------------------------------------------------------------------------
function s = pstr(p)
% PSTR  Compact one-line parameter description.  INPUTS: p.  OUTPUTS: s.
f = setdiff(fieldnames(p), {'muStar','contTol'}, 'stable');
s = '';
for k = 1:numel(f), s = [s sprintf('%s=%.4g ', f{k}, p.(f{k}))]; end %#ok<AGROW>
end

function s = admstr(t)
% ADMSTR  yes/NO marker.  INPUTS: t logical.  OUTPUTS: s [char].
if t, s = 'yes'; else, s = 'NO'; end
end
