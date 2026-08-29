function cat_ = build_gto_catalog()
%% Purpose:
%
%   Packages the completed GTO -> TULIP sheets into the shareable compact
%   catalog costate_catalog_gto_tulip.mat, via the family-agnostic packager
%   costate_common/build_costate_catalog_family. Run from MATLAB after the
%   catalog run completes (16/16 sheets, 2,625 accepted entries as of
%   2026-08-29).
%
%   ORIENTATION-AXIS CONVENTION: the GTO departure is a pseudo-family (an
%   algebraic, unpropagated Kepler locus, not a periodic CR3BP orbit --
%   see get_family_orbit.m's 'gto' case). Its family parameter is
%   orientDeg: the argument-of-perigee angle FROM the Earth->Moon line
%   (the synodic x-axis at epoch) TO the perigee direction, positive in
%   the direction of orbital motion (the same prograde/counter-clockwise
%   sense as the CR3BP synodic frame). sma_km and ecc are held fixed
%   across the whole fleet (350 x 35786 km altitude GTO); only orientDeg
%   varies sheet to sheet.
%
%   SHEET KEY RULE (all-reviewer critical fix): because 'gto' has no
%   periodic period, sheets(k).tauDRO / .tau_dep hold orientDeg (DEGREES,
%   not an ND period) -- run_gto_catalog.m stamped this into every sheet's
%   meta.tauDRO at generation time. This function does NOT touch that; it
%   only verifies it survives packaging unchanged (one assert per sheet
%   below). The true GTO Kepler period rides in dep_params/derive only,
%   never in the key.
%
%   DEP_PARAMS COMPLETION (this function's one substantive step): the
%   sheets on disk only recorded meta.depParams.orientDeg (sma_km/ecc were
%   left to get_family_orbit's defaults at solve time, since the whole
%   fleet uses one fixed GTO). For the packaged catalog's dep_params to be
%   a SELF-CONTAINED, sheet-local reconstruction recipe -- not one that
%   silently depends on get_family_orbit's defaults staying put -- this
%   function completes each sheet's depParams to the full {sma_km, ecc,
%   orientDeg} triple (values pulled from get_family_orbit itself, so they
%   are exactly what reconstruction will use) before handing the sheets to
%   the packager. This patches the ON-DISK sheet .mat files in place
%   (meta only; idempotent -- a sheet already carrying sma_km/ecc is left
%   alone, and the zero-patch case WRITES NOTHING) rather than the
%   packager, since the packager's contract is to copy dep_params verbatim
%   from sheet meta for every family. The rewrite itself is ATOMIC (fix
%   round 2, 2026-08-29 review): a multi-day sheet-fleet campaign has no
%   cheap recovery, so a patched sheet is written to a same-directory temp
%   file, verified by loading it back, THEN swapped in via movefile
%   (atomic on the same filesystem) -- a crash/kill/disk-full mid-write
%   now leaves the ORIGINAL sheet untouched instead of a truncated one.
%
%% Inputs:
%
%   none
%
%% Outputs:
%
%  cat_                     struct                  The catalog, also saved
%                                                   to results/
%                                                   costate_catalog_gto_tulip.mat
%
%% Revision History:
%  M. Casey                                                   (c) 08/29/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(here, 'setup_paths.m'));

catDir = fullfile(here, 'results', 'catalog');
sheetGlob = 'gto_o*_Np*.mat';

complete_dep_params(catDir, sheetGlob);

spec = struct( ...
  'glob', sheetGlob, ...
  'name', 'costate_catalog_gto_tulip', ...
  'description', ['Minimum-time low-thrust GTO->Tulip costate CATALOG: ', ...
    'one (GTO orientation x tulip petal count) phasing sheet, thrust ', ...
    'rungs walked 15->5 N. COMPACT format: canonical ND quantities only; ', ...
    'derived values reconstructed via .derive. Each entry z8 = ', ...
    '[lambda_r(3); lambda_v(3); lambda_m; tf] is a converged PMP solution ', ...
    'accepted by pumpkyn.cr3bp.tfMin. Departure is a GTO PSEUDO-family ', ...
    '(an algebraic, unpropagated Kepler locus -- see get_family_orbit.m ', ...
    'and dep_params below), fixed at 350 x 35786 km altitude; orientDeg ', ...
    'is the argument-of-perigee angle from the Earth->Moon line (the ', ...
    'synodic x-axis at epoch) to the perigee direction, positive in the ', ...
    'direction of orbital motion (CR3BP synodic prograde sense). NOTE ', ...
    '(review fix): because GTO has no orbital PERIOD to key sheets by, ', ...
    'sheets(k).tauDRO/.tau_dep hold orientDeg IN DEGREES here, not an ND ', ...
    'period -- the only catalog in this family where that legacy field ', ...
    'is not a period. The true Kepler period lives in dep_params/derive.'], ...
  'provenance', ['GTO catalog run (run_gto_catalog): warm-recipe fleet, ', ...
    '2026-08-27..29, all 16 sheets (4 orientations {0,90,180,270} deg x ', ...
    '4 petal counts {5,7,9,12}); two-point top-rung multistart per cell ', ...
    '(tf0=0.30 warm pass then cold-tf0 mop-up, first-success-wins); three ', ...
    'gates per entry (ms residual, flown arrival, tfMin acceptance). ', ...
    'Arrival tulips pm = -1 branch; the pm = +1 branch is the exact ', ...
    'z-mirror (CR3BP symmetry). Casey/Koblick 2026-08.'], ...
  'depReconstruction', ['gto: [tau,rv] = get_family_orbit(''gto'', ', ...
    'dep_params) -> algebraic Kepler locus at fixed orientation (NO ', ...
    'propagation); state at phase fraction f: interp1(tau, rv, ', ...
    'mod(f,1)*tau(end)); tulip: getTulip(2*pi*(Np-2)/(Np-1), Np, ', ...
    'pm)->cont_np->prop; state at phase f: interp1(t, rv, ', ...
    'mod(f,1)*t(end), ''spline'')']);

outMat = fullfile(here, 'results', 'costate_catalog_gto_tulip.mat');
cat_ = build_costate_catalog_family(catDir, outMat, spec);
end

% ------------------------------------------------------------------------
function complete_dep_params(catDir, sheetGlob)
% COMPLETE_DEP_PARAMS  Patch each on-disk GTO sheet's meta.depParams to the
% full {sma_km, ecc, orientDeg} triple, idempotently and ATOMICALLY: an
% incomplete sheet is written to a same-directory temp file, verified by
% loading it back, then swapped in via movefile (atomic on the same
% filesystem) -- never a load->mutate->save onto the original path. A
% sheet already complete is skipped entirely (no temp file, no write).
% INPUTS: catDir char (sheet folder); sheetGlob char (file pattern).
% OUTPUTS: none (sheets rewritten on disk, atomically, when incomplete).
files = dir(fullfile(catDir, sheetGlob));
files = files(~contains({files.name}, 'pre_'));
assert(~isempty(files), 'complete_dep_params: no sheet files matching %s in %s', ...
    sheetGlob, catDir);
for kf = 1:numel(files)
    fpath = fullfile(files(kf).folder, files(kf).name);
    S = load(fpath);
    dp = S.meta.depParams;
    if isfield(dp, 'sma_km') && isfield(dp, 'ecc')
        continue                       % already complete -- idempotent
    end
    assert(strcmpi(S.meta.depFamily, 'gto'), ...
        'complete_dep_params:family', ...
        'sheet %s: depFamily is ''%s'', expected ''gto''', ...
        files(kf).name, S.meta.depFamily);
    [~, ~, info] = get_family_orbit('gto', struct('orientDeg', dp.orientDeg));
    S.meta.depParams = struct('sma_km', info.sma_km, 'ecc', info.ecc, ...
        'orientDeg', info.orientDeg);
    assert(abs(info.orientDeg - S.meta.tauDRO) < 1e-12, ...
        'complete_dep_params:key', ...
        'sheet %s: orientDeg (%.6g) ~= tauDRO key (%.6g)', ...
        files(kf).name, info.orientDeg, S.meta.tauDRO);
    % ATOMIC WRITE (fix round 2, 2026-08-29 review): save to a temp file
    % in the SAME directory (movefile is atomic only same-filesystem),
    % verify the temp loads cleanly and carries the patched field, THEN
    % replace the sheet -- a crash/kill/disk-full mid-save now leaves the
    % ORIGINAL sheet untouched instead of a truncated one:
    tmpDir = fileparts(fpath);
    tmpPath = [tempname(tmpDir) '.mat'];   % unique name, SAME dir as fpath
    save(tmpPath, '-struct', 'S');
    try
        Schk = load(tmpPath, 'meta');
        okTmp = isfield(Schk.meta.depParams, 'sma_km') && ...
                isfield(Schk.meta.depParams, 'ecc') && ...
                abs(Schk.meta.depParams.orientDeg - info.orientDeg) < 1e-12;
    catch
        okTmp = false;
    end
    if ~okTmp
        if isfile(tmpPath), delete(tmpPath); end   % no orphaned temp file
        error('complete_dep_params:tmpVerify', ...
            'sheet %s: temp file failed post-write verification -- original left untouched', ...
            files(kf).name);
    end
    movefile(tmpPath, fpath);              % atomic on the same filesystem
    fprintf('complete_dep_params: patched %s -> sma_km=%.6g ecc=%.6g orientDeg=%g\n', ...
        files(kf).name, info.sma_km, info.ecc, info.orientDeg);
end
end
