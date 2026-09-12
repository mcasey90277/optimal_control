function seed = seed_from_entry(entry, rv0, opts)
%% Purpose:
%
%   ONE library entry becomes ONE multiple-shooting seed. A library entry
%   comes in two shapes and each needs a different construction:
%
%     FILE-BACKED  it already carries junction states (an anchor or sweep
%                  result): reuse them, append the final column the solver
%                  wants, and pin column 1 to the ACTUAL departure state and
%                  the entry's costates -- the stored trajectory began at
%                  whatever endpoint IT was solved from;
%     CATALOG      it carries z8 and no junction states: fly it once and cut
%                  the flight (seed_from_z8), which puts the seed AT a root,
%                  so the polish takes one or two Newton steps.
%
%   This construction existed three times -- inline in transfer_study, as
%   run_dro_tulip's private seed_of, and inside build_arrival_sheet's
%   seeding loop -- and all three now call this.
%
%% Inputs:
%
%  entry                    struct                  Library entry: .z [8x1]
%                                                   ([lambda0(7); tf]) and
%                                                   .Y [14 x K] junction
%                                                   START states, or empty
%                                                   .Y for a catalog entry
%  rv0                      [6 x 1] or [1 x 6]      Departure state the seed
%                                                   must start from
%  opts                     struct (optional)
%   .K [24] segments, used for a CATALOG entry (a file-backed entry brings
%   its own), .Tnd .cnd .muStar -- required to FLY a catalog entry
%
%% Outputs:
%
%  seed                     struct                  .tf .tGrid [1 x K+1]
%                                                   .Y [14 x K+1], the
%                                                   ms_bvp / ms_tfmin seed
%                                                   contract
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, opts = struct(); end
            rv0 = rv0(:);

if ~(isstruct(entry) && isfield(entry, 'z') && numel(entry.z) == 8)
    error('seed_from_entry:z8', 'the entry must carry z = [lambda0(7); tf] (8 elements)');
end
if numel(rv0) < 6
    error('seed_from_entry:rv0', 'the departure state needs six components, got %d', numel(rv0));
end
              z = entry.z(:);
          hasYj = isfield(entry, 'Y') && ~isempty(entry.Y);

if ~hasYj
    %% CATALOG entry: fly z8 once and cut. This needs the physics, and a
    %% caller who did not supply it is asking for a seed that cannot exist.
    need = {'Tnd', 'cnd', 'muStar'};
    for k = 1:numel(need)
        if ~isfield(opts, need{k}) || isempty(opts.(need{k}))
            error('seed_from_entry:physics', ...
                  ['this entry carries z8 and no junction states, so it must be FLOWN: ' ...
                   'opts.%s is required'], need{k});
        end
    end
              K = fieldd(opts, 'K', 24);
           seed = seed_from_z8(z, rv0(1:6), K, opts.Tnd, opts.cnd, opts.muStar);
    return
end

%% FILE-BACKED entry: its own junction states, with the endpoint pinned.
              K = size(entry.Y, 2);
           seed = struct('tf', z(8), 'tGrid', linspace(0, z(8), K+1), ...
                         'Y', [entry.Y, entry.Y(:,end)]);
 seed.Y(1:7,1)  = [rv0(1:6); 1];
 seed.Y(8:14,1) = z(1:7);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, d)
% FIELDD  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
