function term = geo_terminal(mode, par, Lf)
% GEO_TERMINAL  Terminal-condition builder for the GEO target.
%
% mode 'fixed'    - full 6-state rendezvous at longitude Lf on GEO (M0/M1).
% mode 'manifold' - free-longitude insertion manifold (M2+): 5 residuals
%                   [r_z; v_z; ||r||^2-a^2; ||v||^2-mu/a; r.v] = 0. NB the set
%                   also admits the retrograde orbit; the prograde seed selects
%                   the branch (process/DESIGN.md sec 2).
%
% INPUTS:  mode - 'fixed' | 'manifold';  par - kepler_lt_params struct;
%          Lf - GEO longitude [rad] ('fixed' only; pass [] for 'manifold')
% OUTPUTS: term - struct (.type; .rvf/.Lf for fixed; .aGeo/.resid for manifold)
%
% REFERENCES:
%   [1] process/DESIGN.md sec 2 (boundary conditions).
switch lower(mode)
    case 'fixed'
        [rf, vf] = elements_to_cart(1, 0, 0, 0, 0, Lf, par.mu);
        term = struct('type','fixed', 'rvf', [rf; vf], 'Lf', Lf);
    case 'manifold'
        a = 1;
        res = @(rv) [rv(3); rv(6); ...
                     rv(1)^2+rv(2)^2+rv(3)^2 - a^2; ...
                     rv(4)^2+rv(5)^2+rv(6)^2 - par.mu/a; ...
                     rv(1)*rv(4)+rv(2)*rv(5)+rv(3)*rv(6)];
        term = struct('type','manifold', 'aGeo', a, 'resid', res);
    otherwise
        error('geo_terminal:mode', 'unknown mode %s', mode);
end
end
