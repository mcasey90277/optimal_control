function K = status_key()
%% Purpose:
%
%   The legend of the optimality STATUS carried by every library entry and
%   alternative: what each code means and the order of the checks that
%   define it. Stored in the catalog as .status_key so the file explains
%   itself to a recipient.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  K                        struct                  .codes [1x6] .names
%                                                   {1x6} .meaning {1x6}
%                                                   .order .note
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

K = struct();
K.codes = [4 3 2 1 0 -1];
K.names = {'sufficient', 'necessary only', 'conjugate point found', 'neither', 'empty', 'below floor'};
K.meaning = { ...
    'first- and second-order sufficient conditions hold (the full certify_root stack passed)', ...
    'first-order necessary conditions hold (pointwise Pontryagin checks + tfMin witness); sufficiency not established', ...
    'an extremal with a conjugate point on the arc, under the hypotheses (gates + H6) that make it a refutation: NOT locally optimal', ...
    'flies to the target, but the necessary conditions are not established', ...
    'no transfer recorded in this cell', ...
    'not a transfer of the library (no converged root, or it does not fly to the target); never stored'};
K.order = ['certify_root: 1 polish, 2 flight (100 km / 10 m/s), 3 pointwise Pontryagin, 4 tfMin witness, ' ...
           '5 coarse conjugate test, 6 hypothesis gates (min|lam_v|, Q_mt, dim S, X2, lift margin), 7 H6, 8 dense conjugate scan'];
K.note = ['A status records what was ESTABLISHED, not what might be true: a check that could not run gives the ' ...
          'highest tier actually established, with the reason. An orbit transfer does not have to be optimal to be useful.'];
end
