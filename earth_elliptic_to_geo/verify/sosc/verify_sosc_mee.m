function sosc = verify_sosc_mee(saved_or_path, opts)
% VERIFY_SOSC_MEE  NLP-level SOSC local-minimum certificate for a saved
% certified MEE min-fuel row. Orchestrates recover -> KKT re-check -> active
% set -> inertia -> verdict + tiered-gate status.
%
% INPUTS:
%   saved_or_path - a sosc_load_row struct OR a .mat path [char]
%   opts - optional struct: .tol (override sosc_defaults)
% OUTPUTS:
%   sosc - struct per process/DESIGN_sosc.md sec 5 (.verdict .reason .status
%          .drift .sign .kkt .active .inertia .redMinEig .thresholds .meta)
% REFERENCES: process/DESIGN_sosc.md secs 4-5.
if nargin<2, opts=struct(); end
tol = optdef(opts,'tol',sosc_defaults());
if ischar(saved_or_path)||isstring(saved_or_path), saved = sosc_load_row(char(saved_or_path));
else, saved = saved_or_path; end

R = sosc_recover_kkt(saved, tol);
sosc = struct('thresholds',tol,'drift',NaN,'sign',NaN, ...
    'kkt',[],'active',[],'inertia',[],'red',[],'nFlat',NaN,'redMinEig',NaN, ...
    'meta',struct('thrustN',saved.thrustN,'tag',saved.tag,'when',datestr(now)));
if ~R.recoverOK
    sosc.verdict='ERROR'; sosc.reason=sprintf('warm re-solve failed: %s',R.ipoptStatus);
    sosc.status='certified-feasibility+sosc-inconclusive'; return;
end
sosc.drift = R.drift;
K  = sosc_kkt_residual(R, tol);   sosc.sign=K.sign;  sosc.kkt=K;
AS = sosc_active_set(R, K, tol);  sosc.active=AS;
IN = sosc_inertia(R.H, AS.A, tol);sosc.inertia=IN;
sosc.red=IN.red; sosc.nFlat=IN.red.nzero;
v  = sosc_decide(K, AS, IN);
sosc.verdict=v.verdict; sosc.reason=v.reason; sosc.status=v.status;
sosc.meta.n=R.n; sosc.meta.m=R.m; sosc.meta.m_active=AS.m_active;
if R.drift >= tol.drift
    warning('verify_sosc_mee:drift','warm re-solve drift %.2e >= %.1e (certifying the re-converged point)', R.drift, tol.drift);
end
end
