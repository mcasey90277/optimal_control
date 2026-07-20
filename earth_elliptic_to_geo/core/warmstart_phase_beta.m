function U = warmstart_phase_beta(X, sigmaDst, dLdst, par, thrRow)
% WARMSTART_PHASE_BETA  Phase-correct RTN thrust-direction (beta) for a warm
% start, replacing the sigma-linear interpolation of beta that aliases across
% rungs (external review, Gemini, 2026-07-19).
%
% The optimal beta oscillates with the ORBITAL PHASE (true longitude L mod
% 2*pi), so interpolating a source beta sequence against the transfer fraction
% sigma keeps the SOURCE revolution frequency: onto a finer-rev rung (more revs)
% the interpolated beta is phase-aliased -- e.g. a 0.5 N -> 0.2 N warm chain has
% 2.5x the revolutions, so a sigma-interpolated beta commands apogee-pointing
% thrust at perigee. This function instead recomputes beta from the SAME
% tangential steering law the seed uses (mee_seed.m>local_beta): the unit RTN
% projection of the local velocity direction, evaluated FRESH at each target
% node's own state and true longitude. Phase-correct by construction, no aliasing.
% (Coasts leave beta free in the NLP, so a tangential guess there is harmless.)
%
% INPUTS:
%   X        - warm-start MEE states [P;ex;ey;hx;hy;m;t] at each node [7x(N+1)]
%   sigmaDst - target sigma grid, [0,1] [(N+1)x1]
%   dLdst    - target total true-longitude span [scalar, ND]
%   par      - kepler_lt_params struct (uses par.mu)
%   thrRow   - throttle row to keep on U row 4 [1x(N+1)]
%
% OUTPUTS:
%   U - warm-start controls [beta(3); thr] [4x(N+1)] with rows 1-3 the
%       phase-correct tangential beta and row 4 = thrRow
%
% REFERENCES:
%   [1] mee_seed.m>local_beta (the tangential steering law reused here).
%   [2] external review (Gemini 3.1 Pro), 2026-07-19: beta sigma-interp aliasing.
N1 = size(X, 2);
sg = sigmaDst(:);
U  = zeros(4, N1);
for k = 1:N1
    Lk = pi + sg(k)*dLdst;                        % this node's true longitude
    [r, v] = elements_to_cart(X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), Lk, par.mu);
    [rhat, that, ~] = rtn_frame(r, v);
    vhat = v / norm(v);
    b = [dot(vhat, rhat); dot(vhat, that); 0];    % v-hat in RTN (normal comp 0)
    U(1:3,k) = b / norm(b);
    U(4,k)   = thrRow(k);
end
end
