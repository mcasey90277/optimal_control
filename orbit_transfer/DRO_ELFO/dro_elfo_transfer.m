% dro_elfo_transfer
%   One minimum-time low-thrust transfer from a distant retrograde orbit (DRO) to an elliptical lunar
%   frozen orbit (ELFO) in the Earth-Moon CR3BP: solve it (the thrust ladder walks down from 15 N to
%   thrust_N), re-fly it independently, check every gate, print one verdict table, save the result and the
%   library-format entry, and draw the figure. Edit the knobs below and run. The target orbit is one
%   struct handed to pair_transfer, so any family it knows can replace the ELFO (see the commented tulip
%   line). A struct dro_elfo_opts in the workspace overrides knobs field by field for ONE run (test seam):
%   the script clears it when it has read it, so a later plain run uses the knobs as written here.

%% knobs
dro_tau      = 1;                     % DRO period [ND], 1 = 4.43 d
elfo         = struct('source', 'im_elfo_optimum', 'sat', 1);   % the IM ELFO satellite 1 (shape, plane, phase)
% elfo       = struct('source', 'elements', 'oe', [12000 0.69 56.5 90]);   % your own [sma km, ecc, inc deg, argp deg]
raan_deg     = NaN;                   % ELFO node angle in the ROTATING frame at the arrival instant [deg]
                                      % (NaN with an IM satellite = that plane's RAAN value, used as this angle:
                                      % a choice of plane orientation, not an epoch; NaN with elements = 0)
sD           = 0;                     % departure phase on the DRO [fraction of its period]
sA           = NaN;                   % arrival phase [fraction] (NaN with an IM satellite = its mean anomaly / 360; otherwise 0)
target       = [];                    % [] = build the ELFO from elfo, raan_deg (ignored when target is given)
% target     = struct('family', 'tulip', 'params', struct('Np', 7, 'pm', -1));   % any get_family_orbit family; then sA is the tulip phase
thrust_N     = 5;                     % thrust [N]
m0_kg        = 150;                   % initial mass [kg]
isp_s        = 1710;                  % specific impulse [s]
floor_km     = 500;                   % minimum altitude above the Moon [km]
max_iter     = 6000;                  % IPOPT iterations per rung (1 = force a failed solve)
make_plot    = true;                  % draw the figure (and save the PNG next to the .mat)
make_movie   = false;                 % movie (not built yet)
out_dir      = '';                    % '' = results/ beside this script
solver       = @pair_transfer;        % test seam, leave as is

knob_names = {'dro_tau', 'elfo', 'raan_deg', 'sD', 'sA', 'target', 'thrust_N', 'm0_kg', 'isp_s', ...
              'floor_km', 'max_iter', 'make_plot', 'make_movie', 'out_dir', 'solver'};
if exist('dro_elfo_opts', 'var')
    for ff = fieldnames(dro_elfo_opts).'
        if ~ismember(ff{1}, knob_names), error('dro_elfo_transfer:knob', 'unknown dro_elfo_opts field ''%s''', ff{1}); end
        eval([ff{1} ' = dro_elfo_opts.(ff{1});']);
    end
    clear dro_elfo_opts ff
end
here = fileparts(mfilename('fullpath'));
if isempty(out_dir), out_dir = fullfile(here, 'results'); end
setup_dro_elfo_paths();

%% build the two orbits
dep = struct('family', 'dro', 'params', struct('tau', dro_tau));
if isempty(target)
    if ~(isstruct(elfo) && isfield(elfo, 'source') && ismember(elfo.source, {'im_elfo_optimum', 'elements'}))
        error('dro_elfo_transfer:elfo', 'elfo needs .source = ''im_elfo_optimum'' (with .sat) or ''elements'' (with .oe)');
    end
    if strcmp(elfo.source, 'im_elfo_optimum')
        [raan_opt, m0_opt, plane_of, oe] = im_elfo_optimum();
        if ~(isfield(elfo, 'sat') && isscalar(elfo.sat) && any(elfo.sat == 1:numel(plane_of)))
            error('dro_elfo_transfer:elfo', 'elfo.sat must be an IM satellite number 1 to %d', numel(plane_of));
        end
        if isnan(raan_deg), raan_deg = raan_opt(plane_of(elfo.sat)); end
        if isnan(sA),       sA       = m0_opt(elfo.sat) / 360;       end
    else
        if ~(isfield(elfo, 'oe') && numel(elfo.oe) == 4)
            error('dro_elfo_transfer:elfo', 'elfo.oe must be [sma_km ecc inc_deg argp_deg]');
        end
        oe = elfo.oe;
        if isnan(raan_deg), raan_deg = 0; end
    end
    target = struct('family', 'elfo', 'params', struct('oe', oe, 'raanDeg', raan_deg));
end
if isnan(sA), sA = 0; end
fprintf('transfer: DRO (period %g ND, phase %.4f) -> %s (phase %.4f), %g N, %g kg, Isp %g s\n', ...
        dro_tau, sD, target.family, sA, thrust_N, m0_kg, isp_s);

%% solve
% the ladder walks 15 N -> ... -> thrust_N (minutes); every gate is applied to the independent re-flight
T = solver(dep, target, sD, sA, thrust_N, struct('m0kg', m0_kg, 'ispS', isp_s, ...
           'floorKm', floor_km, 'maxIter', max_iter));
% T.ok, T.reason           solved, and why not
% T.tf_nd, T.tf_d          minimum time [ND, d]
% T.mf_kg, T.dV_kms        final mass [kg], delta-V [km/s]
% T.flight                 .t .Y (r v mass-fraction lambda) .alpha along the re-flown transfer
% T.gates                  measured verdicts, .pass, .failed
% T.entry                  the transfer in library-entry format
if ~T.ok
    fprintf('NOT SOLVED: %s\n', T.reason);
    return
end

%% verdicts
fprintf('%s\n', strjoin(plot_pair_transfer('verdict', T), newline));

%% save
if ~isfolder(out_dir), mkdir(out_dir); end
mat_file = fullfile(out_dir, 'dro_elfo_transfer.mat');
knobs = struct('dro_tau', dro_tau, 'target', target, 'sD', sD, 'sA', sA, 'thrust_N', thrust_N, ...
               'm0_kg', m0_kg, 'isp_s', isp_s, 'floor_km', floor_km, 'max_iter', max_iter);
save(mat_file, 'T', 'knobs');
% captured entries: a cell array, one library-format entry per cell (entries may differ in fields)
cap_file = fullfile(out_dir, 'dro_elfo_captured_entries.mat');
entries = {};
if isfile(cap_file)
    load(cap_file, 'entries');
    if isstruct(entries), entries = num2cell(entries(:).'); end
end
entries{end + 1} = T.entry;
save(cap_file, 'entries');
fprintf('saved: %s\ncaptured entries: %d in %s\n', mat_file, numel(entries), cap_file);

%% plot
if make_movie, fprintf('movie: not built yet\n'); end
if make_plot
    png_file = fullfile(out_dir, 'dro_elfo_transfer.png');
    plot_pair_transfer(T, struct('file', png_file, 'title', sprintf('DRO to %s, %g N: t_f %.2f d, m_f %.1f kg', ...
                       upper(target.family), thrust_N, T.tf_d, T.mf_kg)));
    fprintf('figure: %s\n', png_file);
end
