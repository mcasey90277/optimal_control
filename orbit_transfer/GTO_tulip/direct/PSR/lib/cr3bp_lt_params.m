function p = cr3bp_lt_params(thrust_N, m0_kg, Isp_s)
% CR3BP_LT_PARAMS  Earth-Moon CR3BP + low-thrust nondimensional parameters.
%
% Centralizes the constants and derived low-thrust quantities that were
% previously copy-pasted into every driver. All ND (nondimensional) values use
% the standard Earth-Moon CR3BP scale set (muStar, lStar, tStar).
%
% INPUTS:
%   thrust_N - max thrust [N, scalar]     (default 0.025 = 25 mN)
%   m0_kg    - initial wet mass [kg]       (default 15)
%   Isp_s    - specific impulse [s]        (default 2100)
%
% OUTPUTS:
%   p - struct with fields:
%       muStar,lStar,tStar - Earth-Moon CR3BP scale set [scalars]
%       m0kg,Isp           - copies of the mass / Isp inputs
%       g0                 - standard gravity, ND [scalar]
%       c                  - exhaust velocity Isp*g0, ND [scalar]
%       Tmax               - max thrust acceleration at m0, ND [scalar]
%
% REFERENCES:
%   [1] Standard Earth-Moon CR3BP normalization (e.g. Koon, Lo, Marsden, Ross).

if nargin < 1 || isempty(thrust_N), thrust_N = 0.025; end
if nargin < 2 || isempty(m0_kg),    m0_kg    = 15;    end
if nargin < 3 || isempty(Isp_s),    Isp_s    = 2100;  end

p.muStar = 0.012150585609624;
p.lStar  = 389703.264829278;      % km
p.tStar  = 382981.289129055;      % s
p.m0kg   = m0_kg;
p.Isp    = Isp_s;
p.g0     = 9.80665 * p.tStar^2 / (1000*p.lStar);        % ND
p.c      = (Isp_s / p.tStar) * p.g0;                    % ND exhaust velocity
p.Tmax   = (thrust_N / m0_kg) * p.tStar^2 / (p.lStar*1000);  % ND accel
end
