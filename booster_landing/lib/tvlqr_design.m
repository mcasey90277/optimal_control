function ctrl = tvlqr_design(sol, P, opts)
% TVLQR_DESIGN  Time-varying LQR about a PDG guidance trajectory.
%
% Linearizes pdg_dynamics along (x*(t), T*(t)) and integrates the
% differential Riccati equation backward:
%   -Pdot = A'P + PA - PBR^{-1}B'P + Q,   P(tf) = Qf
% Gains K(t) = R^{-1} B(t)' P(t), stored on a dense grid; controller is
%   T_cmd(t) = T*(t) - K(t) (x - x*(t))  (saturation applied by the sim).
%
% INPUTS:
%   sol  - collocation solution (Task 3 interface)
%   P    - booster_params
%   opts - (optional) .Q .R .Qf weight overrides (7x7, 3x3, 7x7)
% OUTPUTS:
%   ctrl - .tgrid(1xM) .K(3x7xM) .Pt(7x7xM) .xnom(@t->7x1) .Tnom(@t->3x1)
%
% REFERENCES:
%   [1] Anderson & Moore, "Optimal Control: Linear Quadratic Methods."
if nargin < 3, opts = struct(); end
if ~isfield(opts,'Q'),  opts.Q  = diag([1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 0]); end
if ~isfield(opts,'R'),  opts.R  = 1e-10*eye(3);                            end
if ~isfield(opts,'Qf'), opts.Qf = diag([1e-2 1e-2 1e-2 1 1 1 0]);          end

%% Nominal interpolants (nodes+midpoints for thrust, nodes for state):
Nn = size(sol.X,2) - 1;  h = sol.tf/Nn;
tU = sort([sol.t, sol.t(1:end-1) + h/2]);
TU = zeros(3, numel(tU));  TU(:,1:2:end) = sol.U;  TU(:,2:2:end) = sol.Um;
ctrl.xnom = @(t) interp1(sol.t.', sol.X.', clampt(t, sol.tf), 'pchip').';
ctrl.Tnom = @(t) interp1(tU.',   TU.',    clampt(t, sol.tf), 'pchip').';

%% Backward Riccati on vec(P), dense output grid:
M     = 4*(Nn+1);
ctrl.tgrid = linspace(0, sol.tf, M);
Rinv  = inv(opts.R);
ric   = @(t, pv) ricrhs(t, pv, ctrl, P, opts.Q, Rinv);
oo    = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, PV] = ode45(ric, fliplr(ctrl.tgrid), opts.Qf(:), oo);   % integrates tf->0
PV    = flipud(PV);                                          % re-order 0->tf
ctrl.Pt = zeros(7,7,M);  ctrl.K = zeros(3,7,M);
for k = 1:M
    Pk = reshape(PV(k,:), 7, 7);  Pk = (Pk + Pk.')/2;
    [~, ~, Bk] = pdg_dynamics(ctrl.xnom(ctrl.tgrid(k)), ...
                              ctrl.Tnom(ctrl.tgrid(k)), P);
    ctrl.Pt(:,:,k) = Pk;
    ctrl.K(:,:,k)  = Rinv * (Bk.' * Pk);
end
end

function t = clampt(t, tf), t = min(max(t, 0), tf); end

function pdot = ricrhs(t, pv, ctrl, P, Q, Rinv)
Pk = reshape(pv, 7, 7);
[~, A, B] = pdg_dynamics(ctrl.xnom(t), ctrl.Tnom(t), P);
Pd   = -(A.'*Pk + Pk*A - Pk*B*Rinv*B.'*Pk + Q);
pdot = Pd(:);
end
