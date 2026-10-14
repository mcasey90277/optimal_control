function lines = print_transfer_summary(S, opts)
%% Purpose:
%
%   The transfer numbers every consumer prints, printed one way: time of
%   flight, Delta-V, propellant, and -- when the caller has the flight
%   itself -- how many samples it took and how close it came to each
%   primary.
%
%   Serves BOTH the study script (which holds a fly_transfer flight) and
%   run_dro_tulip (which holds a certify_root certificate), because they
%   were printing the same three numbers in two formats. Each caller keeps
%   its own extra lines: the script's lambda_0, the front door's gates,
%   witness and source.
%
%  ASSUMPTIONS / NOTES:
%
% • It prints what the struct HAS. A certificate carries t_f in days only,
%   a flight carries it in ND too; the flight line appears only when the
%   admissibility verdict and its clearances are there.
%
%% Inputs:
%
%  S                        struct                  A fly_transfer flight or
%                                                   a certify_root result:
%                                                   needs .dvKms
%                                                   .propellantKg and one of
%                                                   .tfDays / .tfNd;
%                                                   optional .nSamples
%                                                   .admissibility (.reason
%                                                   .moonKm .earthKm)
%  opts                     struct (optional)
%   .prefix ['   '] line prefix, .quiet [false] build the lines without
%   printing them
%
%% Outputs:
%
%  lines                    cellstr                 the lines, printed
%                                                   unless quiet
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
         prefix = fieldd(opts, 'prefix', '   ');
          quiet = fieldd(opts, 'quiet', false);
         tStar_ = fieldd(S, 'tStar', 382981.289129055);

for f = {'dvKms', 'propellantKg'}
    if ~isfield(S, f{1}) || isempty(S.(f{1}))
        error('print_transfer_summary:fields', 'the struct carries no %s', f{1});
    end
end
if isfield(S, 'tfDays') && ~isempty(S.tfDays)
         tfDays = S.tfDays;
elseif isfield(S, 'tfNd') && ~isempty(S.tfNd)
         tfDays = S.tfNd*tStar_/86400;
else
    error('print_transfer_summary:fields', 'the struct carries no time of flight');
end

          lines = {};
if isfield(S, 'tfNd') && ~isempty(S.tfNd)
    lines{end+1} = sprintf('%st_f = %.6f ND (%.4f d)   Delta-V = %.4f km/s   propellant %.2f kg', ...
                           prefix, S.tfNd, tfDays, S.dvKms, S.propellantKg);
else
    lines{end+1} = sprintf('%st_f = %.4f d   Delta-V = %.4f km/s   propellant %.2f kg', ...
                           prefix, tfDays, S.dvKms, S.propellantKg);
end
if isfield(S, 'admissibility') && isstruct(S.admissibility) && isfield(S.admissibility, 'moonKm')
    lines{end+1} = sprintf('%sflight: %d samples, %s; closest approach %.0f km (Moon) / %.0f km (Earth)', ...
                           prefix, fieldd(S, 'nSamples', NaN), S.admissibility.reason, ...
                           S.admissibility.moonKm, S.admissibility.earthKm);
end

if ~quiet
    for k = 1:numel(lines), fprintf('%s\n', lines{k}); end
end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, d)
% FIELDD  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
