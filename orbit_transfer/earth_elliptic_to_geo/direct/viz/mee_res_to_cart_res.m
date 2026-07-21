function cartRes = mee_res_to_cart_res(Xmee, Umee, dL, sigma, thrustN, ctf, mu, nDense)
% MEE_RES_TO_CART_RES  Convert an MEE/L-domain min-fuel solution into the
% Cartesian (inertial) results layout that transfer_movie.m consumes.
%
% The MEE solver (casadi_lt_mee.m) stores its state in elements
% [P;ex;ey;hx;hy;m;t] with true longitude L = pi + sigma*dL as the
% independent variable, and its control in the local RTN frame
% [beta(3);thr]. transfer_movie.m, however, was written for the Cartesian/
% Sundman solver and expects a 9-row inertial state [r(3);v(3);m;t;cScale]
% and a 4-row inertial control [alpha(3);thr]. This adapter bridges the two
% by (1) reconstructing inertial (r,v) at every node via elements_to_cart at
% L_k = pi + sigma_k*dL, and (2) rotating the RTN thrust direction beta into
% the inertial frame, alpha = R_{RTN->ECI} * beta, using the same per-node
% RTN triad (rhat, that, nhat) that run_transfer_mee.m's own cross-
% formulation reconstruction check builds. The output plugs straight into
% transfer_movie.m (save it as `res` and pass the file path).
%
% INPUTS:
%   Xmee    - MEE states [P;ex;ey;hx;hy;m;t] at each node [7 x (N+1)]
%   Umee    - MEE controls [beta_RTN(3);thr] at each node [4 x (N+1)]
%   dL      - total true-longitude span L(end)-pi [scalar]
%   sigma   - uniform node parameter 0->1 [(N+1) x 1]
%   thrustN - max thrust [N] (echoed into cartRes.cfg for the movie title)
%   ctf     - t_f / t_f,min ratio [scalar] (echoed into cartRes.cfg)
%   mu      - gravitational parameter in the solution's units [scalar, =1 ND]
%   nDense  - OPTIONAL render-densification factor [scalar int, default 1].
%             nDense=1 reproduces the original node-only output byte-for-byte.
%             nDense>1 inserts nDense-1 intermediate points per original segment
%             by shape-preserving (pchip) interpolation of the MEE elements onto
%             a fine sigma grid and re-mapping to Cartesian at the fine true
%             longitudes L = pi + sigma*dL. Because the equinoctial elements
%             drift SLOWLY while L sweeps 2*pi per revolution, this recovers the
%             smooth orbit that the coarse mesh (~8 nodes/rev at deep rungs)
%             renders as a polygon. Throttle is held piecewise-constant ('previous')
%             so the bang-bang switch structure and switch COUNT are preserved;
%             beta is pchip-interpolated and renormalized (arrow display only).
%
% OUTPUTS:
%   cartRes - struct matching transfer_movie.m's expectations:
%             .cfg.thrustN .cfg.ctf
%             .fuel.X [9 x (N+1)] = [r(3); v(3); m; t; 0]  (row 9 unused)
%             .fuel.U [4 x (N+1)] = [alpha_ECI(3); thr]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/run_transfer_mee.m>check_reconstruction
%       (the per-node RTN triad + elements_to_cart reconstruction reused here).
%   [2] earth_elliptic_to_geo/transfer_movie.m (consumer of this layout).
%   [3] earth_elliptic_to_geo/elements_to_cart.m (algebraic MEE->(r,v) map).

Nn = size(Xmee, 2);
sigma = sigma(:);
assert(numel(sigma) == Nn, 'mee_res_to_cart_res:sizeMismatch', ...
    'sigma has %d entries but Xmee has %d columns', numel(sigma), Nn);
if nargin < 8 || isempty(nDense), nDense = 1; end

% Optional render densification: interpolate the (slowly-varying) elements onto
% a fine sigma grid; L is re-mapped analytically at each fine node so the fast
% orbital sweep is resolved. nDense=1 leaves Xmee/Umee untouched (byte-identical).
if nDense > 1
    sf   = linspace(sigma(1), sigma(end), (Nn-1)*nDense + 1).';
    Xf   = zeros(7, numel(sf));
    for rr = 1:7, Xf(rr,:) = interp1(sigma, Xmee(rr,:), sf, 'pchip'); end
    bt   = interp1(sigma, Umee(1:3,:).', sf, 'pchip').';   % [3 x Nf]
    bt   = bt ./ max(vecnorm(bt,2,1), realmin);            % renormalize beta
    th   = interp1(sigma, Umee(4,:), sf, 'previous').';    % hold throttle (bang-bang)
    Xmee = Xf;  Umee = [bt; th(:).'];  sigma = sf;  Nn = numel(sf);
end

Xc = zeros(9, Nn);
Uc = zeros(4, Nn);
for k = 1:Nn
    Lk = pi + sigma(k)*dL;
    [rk, vk] = elements_to_cart(Xmee(1,k), Xmee(2,k), Xmee(3,k), ...
                                Xmee(4,k), Xmee(5,k), Lk, mu);
    % inertial RTN triad at this node (radial, transverse, normal)
    [rhat, that, nhat] = rtn_frame(rk, vk);
    Rrtn2eci = [rhat, that, nhat];           % columns = R, T, N (inertial)

    beta = Umee(1:3, k);
    alphaEci = Rrtn2eci * beta;              % inertial thrust direction

    Xc(1:3, k) = rk;
    Xc(4:6, k) = vk;
    Xc(7,   k) = Xmee(6, k);                 % mass
    Xc(8,   k) = Xmee(7, k);                 % physical time
    Xc(9,   k) = 0;                          % cScale slot (unused by the movie)
    Uc(1:3, k) = alphaEci;
    Uc(4,   k) = Umee(4, k);                 % throttle
end

cartRes = struct();
cartRes.cfg  = struct('thrustN', thrustN, 'ctf', ctf);
cartRes.fuel = struct('X', Xc, 'U', Uc);
end
