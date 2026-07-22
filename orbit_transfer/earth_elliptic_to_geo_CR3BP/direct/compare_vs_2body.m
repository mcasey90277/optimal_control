% COMPARE_VS_2BODY  Per-rung comparison: CR3BP certified min-fuel solves
% vs. the 2-body certified Table-3 baseline.
%
% Phase-1, no-solve script (runs in seconds). Scans results/ for every
% Task-5 certified artifact (minfuel_cr3bp_T<thrustTag>N_phi<phiTag>.mat),
% and for each one prints/writes a row comparing:
%   - CR3BP m_f_kg (this campaign's certified solve) vs the 2-body certified
%     m_f_kg (table3_certified.m), as Delta_m_f in kg and in percent
%   - switch counts from BOTH sides, explicitly labeled as nodal counts with
%     the P0 mesh-band caveat (neither side is an independently-verified
%     converged switch-mesh count at these rungs -- see
%     table3_certified.m's DEEP-RUNG SWITCH-COUNT CAVEAT and
%     solve_cr3bp_minfuel.m's own switches docstring)
%     NOTE: at 10 N the printed "sw 2b/cr3bp" column is bare integers
%     (e.g. "19/19") with the mesh-band caveat carried only in the footnote
%     below -- acceptable at 10 N because table3_certified.m's
%     under-resolution caveat applies only to the 0.2/0.1 N rungs (counts
%     ARE mesh-adequate here). Spec sec 8 gate 4's literal "never bare
%     integers" wording will need the swStr format itself to carry a band
%     marker (e.g. "~19/~19" or an explicit +/-N) once the deep rungs
%     (1/0.2/0.1 N, TODO.md) are added to this table.
%   - the NLP defect (best.maxDefect) as a certification sanity check
%   - t_f in days and lunar months (t_f = 1.5*tfmin, campaign's c_tf,
%     identical convention on both sides per spec D4)
%   - the sanity_bound.m null-model tide/authority ratio [%], recomputed
%     here (not re-loaded) via lunar_params so this table is self-contained
%
% This is the spec sec-8 gate-4 artifact: "did the certified CR3BP solve(s)
% beat, match, or fall short of the sec-7 null-hypothesis prediction, with
% the switch-count caveat explicit." No solving happens here -- this script
% only loads existing .mat artifacts and the certified 2-body lookup table.
%
% INPUTS:  none (script; globs results/minfuel_cr3bp_*.mat internally)
% OUTPUTS: none (prints an aligned table to stdout; writes
%          orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/results/compare_vs_2body.md)
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md
%       sec 4 (D4: comparison convention), sec 7 (null hypothesis), sec 8
%       (gate 4).
%   [2] orbit_transfer/earth_elliptic_to_geo/direct/reproduce/table3_certified.m
%       (certified 2-body per-rung m_f/switches/revs/tfmin -- source of truth,
%       including its own DEEP-RUNG SWITCH-COUNT CAVEAT).
%   [3] orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/solve_cr3bp_minfuel.m
%       (Task 5; produces the minfuel_cr3bp_*.mat artifacts this script reads).
%   [4] orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/sanity_bound.m
%       (Task 3; the null-model tide/authority ratio, recomputed here per rung
%       actually solved rather than the full certified-rung sweep).

here = fileparts(mfilename('fullpath'));
addpath(here);
setup_paths;

resDir  = fullfile(here, 'results');
pattern = fullfile(resDir, 'minfuel_cr3bp_*.mat');
files   = dir(pattern);

assert(~isempty(files), 'compare_vs_2body:noArtifacts', ...
    'No minfuel_cr3bp_*.mat artifacts found in %s -- run solve_cr3bp_minfuel first', resDir);

c_tf      = 1.5;      % campaign's t_f/t_fMin factor (spec D4; matches sanity_bound.m)
monthDays = 27.32;    % sidereal lunar month [days]

n = numel(files);
rows = cell(n, 1);

for k = 1:n
    S = load(fullfile(files(k).folder, files(k).name));
    best = S.best;
    fp   = S.fp;
    prov = S.provenance;

    cert2body = table3_certified(fp.thrustN);

    par  = kepler_lt_params(fp.thrustN, fp.m0kg, fp.ispS);
    pert = lunar_params(par, fp.phi0, fp.gain);

    row.file        = files(k).name;
    row.thrustN     = fp.thrustN;
    row.phi0        = fp.phi0;
    row.mf_2body    = cert2body.m_f_kg;
    row.mf_cr3bp    = best.m_f_kg;
    row.dmf_kg      = best.m_f_kg - cert2body.m_f_kg;
    row.dmf_pct     = 100 * row.dmf_kg / cert2body.m_f_kg;
    row.sw_2body    = cert2body.switches;
    row.sw_cr3bp    = best.switches;
    row.maxDefect   = best.maxDefect;
    row.epsReached  = best.epsReached;
    row.certified   = best.certified;

    tfDays          = c_tf * cert2body.tfmin * par.TU_s / 86400;
    row.tfDays      = tfDays;
    row.tfMonths    = tfDays / monthDays;

    authority       = fp.thrustN / fp.m0kg;                      % [m/s^2]
    tide            = 2 * pert.muM * 1 / pert.DM^3 * par.AU_ms2; % [m/s^2] at r=1 LU
    row.ratioPct    = 100 * tide / authority;

    row.dmf_kg_predicted_sign = sign(row.dmf_kg);   % + means Moon HELPS at this phi0

    rows{k} = row; %#ok<AGROW>

    % Sanity echo tying dmf back to solve_cr3bp_minfuel's own provenance
    % (must match to machine precision -- both computed the same subtraction).
    if abs(row.dmf_kg - prov.dmf_kg) > 1e-9
        warning('compare_vs_2body:dmfMismatch', ...
            '%s: recomputed dmf_kg=%.6f differs from provenance.dmf_kg=%.6f', ...
            files(k).name, row.dmf_kg, prov.dmf_kg);
    end
end

%% Print aligned table to stdout:
hdr = sprintf('%6s %6s %11s %11s %10s %8s %14s %10s %10s %11s %10s', ...
    'T [N]', 'phi0', 'mf_2body', 'mf_cr3bp', 'dmf [kg]', 'dmf [%]', ...
    'sw 2b/cr3bp*', 'defect', 't_f [d]', 't_f [mo]', 'ratio [%]');
sep = repmat('-', 1, numel(hdr));
fprintf('%s\n%s\n', hdr, sep);
for k = 1:n
    r = rows{k};
    swStr = sprintf('%d/%d', r.sw_2body, r.sw_cr3bp);
    fprintf('%6.3g %6.3g %11.4f %11.4f %10.4f %8.5f %14s %10.3e %10.3f %11.4f %10.4f\n', ...
        r.thrustN, r.phi0, r.mf_2body, r.mf_cr3bp, r.dmf_kg, r.dmf_pct, ...
        swStr, r.maxDefect, r.tfDays, r.tfMonths, r.ratioPct);
end
fprintf(['\n* switch counts are NODAL counts -- mesh-band caveat (P0 protocol): ' ...
    'neither the 2-body nor the CR3BP switch count at these rungs is an\n' ...
    '  independently mesh-converged value (see table3_certified.m''s DEEP-RUNG ' ...
    'SWITCH-COUNT CAVEAT); read as bands, not exact integers.\n\n']);

for k = 1:n
    r = rows{k};
    sgn = 'HELPS';
    if r.dmf_kg < 0, sgn = 'HURTS'; end
    fprintf('%s: certified=%d epsReached=%d maxDefect=%.3e -- Moon %s (dmf=%+.4f kg = %+.5f%%) at phi0=%g rad\n', ...
        r.file, r.certified, r.epsReached, r.maxDefect, sgn, r.dmf_kg, r.dmf_pct, r.phi0);
end
fprintf('\n');

%% Write markdown table:
mdFile = fullfile(resDir, 'compare_vs_2body.md');
fid = fopen(mdFile, 'w');
fprintf(fid, '# CR3BP vs. 2-body: certified min-fuel comparison\n\n');
fprintf(fid, ['Generated by `compare_vs_2body.m` (phase1 T6). Sources: certified\n' ...
    '`table3_certified.m` (2-body) and `minfuel_cr3bp_*.mat` (Task-5 CR3BP\n' ...
    'certified solves) + `lunar_params.m` (null-model ratio, spec sec 7).\n' ...
    'c_tf = %g, lunar month = %g d.\n\n'], c_tf, monthDays);
fprintf(fid, ['| T [N] | phi0 [rad] | m_f 2-body [kg] | m_f CR3BP [kg] | ' ...
    'Delta m_f [kg] | Delta m_f [%%] | switches 2b/CR3BP* | maxDefect | ' ...
    't_f [days] | t_f [lunar months] | tide/authority [%%] |\n']);
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for k = 1:n
    r = rows{k};
    swStr = sprintf('%d/%d', r.sw_2body, r.sw_cr3bp);
    fprintf(fid, '| %g | %g | %.4f | %.4f | %+.4f | %+.5f | %s | %.3e | %.3f | %.4f | %.4f |\n', ...
        r.thrustN, r.phi0, r.mf_2body, r.mf_cr3bp, r.dmf_kg, r.dmf_pct, ...
        swStr, r.maxDefect, r.tfDays, r.tfMonths, r.ratioPct);
end
fprintf(fid, ['\n*Switch counts are NODAL counts with a mesh-band caveat (P0 ' ...
    'protocol) -- neither side is independently mesh-converged at these\n' ...
    'rungs (see `table3_certified.m`''s DEEP-RUNG SWITCH-COUNT CAVEAT); read ' ...
    'as bands, not exact integers.\n\n']);
fprintf(fid, '**Certification status per row:**\n\n');
for k = 1:n
    r = rows{k};
    sgn = 'HELPS';
    if r.dmf_kg < 0, sgn = 'HURTS'; end
    fprintf(fid, ['- `%s`: certified=%d, epsReached=%d, maxDefect=%.3e -- Moon ' ...
        '%s (Delta m_f = %+.4f kg = %+.5f%%) at phi0=%g rad\n'], ...
        r.file, r.certified, r.epsReached, r.maxDefect, sgn, r.dmf_kg, r.dmf_pct, r.phi0);
end
fclose(fid);
fprintf('compare_vs_2body: wrote %s\n', mdFile);
