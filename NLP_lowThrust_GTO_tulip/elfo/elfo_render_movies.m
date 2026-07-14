function elfo_render_movies(sel, mode)
% ELFO_RENDER_MOVIES  Render a control movie (or preview stills) from each saved
% GTO->ELFO min-fuel solution, POST-HOC -- no re-solve. Calls elfo_movie on the
% already-banked minfuel_ELFO_*.mat solution files (e.g. from elfo_batch.sh).
%
% The batch (elfo_run_one) deliberately does not render movies -- it only solves
% and writes the solution .mat. This renders them separately, so the (long)
% solve sweep is never slowed by ~10-15 min/movie rendering.
%
% INPUTS:
%   sel  - 'all' (every results/minfuel_ELFO_tf*_minEps0.mat) | a numeric factor
%          array (e.g. [1.20 1.65]) | [] -> 'all'   [char | 1xM double]
%   mode - 'movie' (MP4+GIF, ~10-15 min each) | 'preview' (3 stills, seconds)
%          [char, default 'movie']
%
% OUTPUTS: none (writes results/movie_ELFO_tf<fTag>_minEps0.{mp4,gif} in 'movie'
%          mode, or ..._{early,mid,late}.png in 'preview' mode, per factor)
%
% REFERENCES:
%   [1] elfo_movie.m (the renderer); [2] gen_elfo_minfuel.m (the solution files);
%   [3] run_elfo_minfuel.m stage 6 (the single-factor movie pattern this reuses).

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
resDir = fullfile(here, 'results');
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
if nargin < 1, sel = []; end
if nargin < 2 || isempty(mode), mode = 'movie'; end

% ELFO one-period backdrop trace (shared by every render)
[~, ~, elfoTrace] = gto_elfo_endpoints(p, struct('point','apolune'));
bg = elfoTrace(:, 1:3);

% --- resolve the solution-file list -----------------------------------------
if isempty(sel) || ((ischar(sel) || isstring(sel)) && strcmpi(char(sel), 'all'))
    d = dir(fullfile(resDir, 'minfuel_ELFO_tf*_minEps0.mat'));
    files = fullfile(resDir, {d.name});
else
    files = {};
    for f = sel(:)'
        fTag = strrep(sprintf('%.3f', f), '.', 'p');
        d = dir(fullfile(resDir, sprintf('minfuel_ELFO_tf%s_sw*_minEps0.mat', fTag)));
        if isempty(d)
            warning('elfo_render_movies:noSol', 'no ELFO min-fuel solution for factor %.3f', f);
            continue
        end
        files{end+1} = fullfile(resDir, d(1).name); %#ok<AGROW>
    end
end
if isempty(files)
    fprintf('elfo_render_movies: no solutions to render\n');  return
end

fprintf('=== ELFO_RENDER_MOVIES (%s): %d solution(s) ===\n', mode, numel(files));
for k = 1:numel(files)
    solFile = files{k};
    try
        L    = load(solFile, 'U', 'factor');
        ss   = L.U(4, :);  nsw = sum(abs(diff(ss > 0.5)));
        fac  = L.factor;
        fTag = strrep(sprintf('%.3f', fac), '.', 'p');
        stem = fullfile(resDir, sprintf('movie_ELFO_tf%s_minEps0', fTag));
        % no switch count in the title -- elfo_movie's panels report it (with its
        % own threshold); a second count here would clash within the same frame.
        ttl  = sprintf('min-fuel GTO\\rightarrowELFO, t_f=%.2fx (bang-bang min-fuel)', fac);
        fprintf('[%d/%d] factor %.3f (%d sw) -> %s  [%s]\n', k, numel(files), fac, nsw, stem, mode);
        elfo_movie(solFile, stem, ttl, mode, bg);
    catch ME
        warning('elfo_render_movies:renderFail', 'render failed for %s: %s', solFile, ME.message);
    end
end
fprintf('ELFO_RENDER_MOVIES DONE\n');
end
