function P = accessibility_plane(T, caps_days)
%% Purpose:
%
%   The (thrust bound, time cap) plane of a sheet table, two ways. RAW
%   COVERAGE: at each rung and cap, how many pairs have an ACCEPTED cell
%   solved at that rung with time of flight under the cap, and how many
%   cells -- what the catalogs hold at that thrust. KNOWN-FEASIBLE
%   ENVELOPE: a trajectory accepted at a lower thrust bound is feasible
%   under a higher one (throttleable thrust), so the envelope at a rung
%   unions every rung at or below its thrust. Thrust here is a capability
%   axis for one spacecraft, not a consumable budget. Rungs are the
%   catalog's, shallow to deep; caps are in days.
%
%% Inputs:
%
%  T                        struct [n x 1]          Sheet table
%                                                   (load_catalog_sheets)
%
%  caps_days                [1 x nC]                Time caps, days,
%                                                   increasing
%
%% Outputs:
%
%  P                        struct                  .rungs_N [1 x nR];
%                                                   .caps_days [1 x nC];
%                                                   .pairs, .cells [nR x
%                                                   nC] raw; .closed {nR
%                                                   x nC} sheet indices
%                                                   in .pairs; .pairs_env,
%                                                   .cells_env [nR x nC]
%                                                   envelope
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the four catalogs on 10 caps between 1 and 4 days
     here = fileparts(fileparts(mfilename('fullpath')));
     addpath(here);
     P = accessibility_plane(load_catalog_sheets([]), linspace(1, 4, 10));
     fprintf('raw pairs (rows = rungs %s N, cols = caps %s d):\n', mat2str(P.rungs_N), mat2str(P.caps_days, 3));
     disp(P.pairs);
     fprintf('envelope pairs:\n');  disp(P.pairs_env);
     return
end

caps  = caps_days(:).';
nC    = numel(caps);
rungs = T(1).rungs_N;
nR    = numel(rungs);
P.rungs_N   = rungs;
P.caps_days = caps;
P.pairs     = zeros(nR, nC);
P.cells     = zeros(nR, nC);
P.closed    = repmat({zeros(1, 0)}, nR, nC);
P.pairs_env = zeros(nR, nC);
P.cells_env = zeros(nR, nC);
for kk = 1:numel(T)
    toDay = T(kk).constants.tStar_s/86400;
    tfA   = T(kk).tf_nd*toDay;
    tfA(~T(kk).accepted) = Inf;                            % unaccepted cells never count
    for rr = 1:nR
        tc = tfA(:, :, rr);  tc = tc(isfinite(tc));        % raw: this rung's accepted cells
        te = min(tfA(:, :, rr:end), [], 3);                % envelope: best over this and every deeper rung, per cell
        te = te(isfinite(te));
        for cc = 1:nC
            nRaw = nnz(tc <= caps(cc));
            nEnv = nnz(te <= caps(cc));
            P.cells(rr, cc)     = P.cells(rr, cc) + nRaw;
            P.cells_env(rr, cc) = P.cells_env(rr, cc) + nEnv;
            if nRaw > 0
                P.pairs(rr, cc)  = P.pairs(rr, cc) + 1;
                P.closed{rr, cc} = [P.closed{rr, cc}, kk];
            end
            if nEnv > 0
                P.pairs_env(rr, cc) = P.pairs_env(rr, cc) + 1;
            end
        end
    end
end
end
