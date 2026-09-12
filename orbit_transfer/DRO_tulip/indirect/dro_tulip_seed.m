function [seed, info] = dro_tulip_seed(sD, sA, op)
%% Purpose:
%
%   WHICH certified solution seeds this phase pair, and is it the right
%   engine? The lookup half of what transfer_study section 4 used to spell
%   out inline: find the pair in the certified library (the anchor and sweep
%   files PLUS the 115-entry 70 mN catalog), verify the whole operating
%   point, and hand back the multiple-shooting seed.
%
%   The construction itself belongs to costate_common/seed_from_entry; this
%   function owns only the lookup and the refusals.
%
%  ASSUMPTIONS / NOTES:
%
% • The phase match is EXACT to 1e-6, never nearest. A neighbouring phase
%   would usually converge, but at one arrival phase this problem carries a
%   ladder of extremals -- one sheet column held 14 candidates, of which the
%   conjugate test refuted 12 -- so a neighbour seed can land on a slower
%   branch and pass every first-order check. Walking there properly is
%   run_dro_tulip's job (continuation, certified at each step).
%
% • The operating point is checked in ALL SIX fields. Two separate
%   identifiers, so "no solution at this phase pair" and "that is a
%   different engine" are distinguishable; the inline version asserted both
%   at once and could not say which had failed.
%
%% Inputs:
%
%  sD, sA                   double                  Departure and arrival
%                                                   phase fractions
%  op                       struct                  The operating point to
%                                                   verify and fly with:
%                                                   .tau .Np .pm .thrustN
%                                                   .ispS .m0kg .rv0
%                                                   (departure state) .Tnd
%                                                   .cnd .muStar, optional
%                                                   .K [24] .tolPhase [1e-6]
%
%% Outputs:
%
%  seed                     struct                  ms seed (.tf .tGrid .Y)
%  info                     struct                  .entry (the library
%                                                   entry) .src .sD .sA
%                                                   .tfDays .nLibrary
%                                                   .fromCatalog
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

           here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(here)), 'costate_common'));
       tolPhase = fieldd(op, 'tolPhase', 1e-6);

%% The operating point the certified library was built at:
          libOp = struct('tau', 1.0, 'Np', 7, 'pm', -1, ...
                         'thrustN', 0.070, 'ispS', 900, 'm0kg', 150);
            tol = struct('tau', 1e-12, 'Np', 0, 'pm', 0, ...
                         'thrustN', 1e-12, 'ispS', 1e-9, 'm0kg', 1e-9);
            fns = fieldnames(libOp);
for k = 1:numel(fns)
              f = fns{k};
    if ~isfield(op, f)
        error('dro_tulip_seed:operatingPoint', 'the operating point is missing %s', f);
    end
    if ~(abs(op.(f) - libOp.(f)) <= tol.(f))
        error('dro_tulip_seed:operatingPoint', ...
              ['this is a different operating point: %s = %g, and the certified library ' ...
               'covers %g (tau = 1 DRO -> 7-petal tulip, pm = -1, 70 mN, Isp 900 s, 150 kg)'], ...
              f, op.(f), libOp.(f));
    end
end

%% The lookup -- exact, never nearest (see the note above):
            lib = dro_tulip_library([], struct('includeCatalog', true));
          match = find(abs([lib.sD] - mod(sD,1)) < tolPhase & ...
                       abs([lib.sA] - mod(sA,1)) < tolPhase, 1);
if isempty(match)
    error('dro_tulip_seed:noSeed', '%s', sprintf( ...
        ['no certified solution at phase pair (%.4f, %.4f). The library holds %d of them\n' ...
         'on the grid sD = k/12, sA = 0.0754 + j/12. For any other pair, walk to it with\n' ...
         '   T = run_dro_tulip(sD, sA, opts)   and study T.'], mod(sD,1), mod(sA,1), numel(lib)));
end

          entry = lib(match);
           seed = seed_from_entry(entry, op.rv0, op);
           info = struct('entry', entry, 'src', entry.src, 'sD', entry.sD, 'sA', entry.sA, ...
                         'tfDays', entry.tfDays, 'nLibrary', numel(lib), ...
                         'fromCatalog', isempty(entry.Y), 'operatingPoint', libOp);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, d)
% FIELDD  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
