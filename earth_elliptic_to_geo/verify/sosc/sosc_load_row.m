function saved = sosc_load_row(matPath)
% SOSC_LOAD_ROW  Normalize a saved certified fuel row (MEE_M2 res-struct or PSR
% out-struct) into the common `saved` struct the SOSC recovery consumes.
%
% INPUTS:
%   matPath - path to a certified .mat (results/MEE_M2_*.mat or
%             results/*_PSR_psr_final.mat) [char]
% OUTPUTS:
%   saved - struct: .sigma[(N+1)x1] .X[7x(N+1)] .U[4x(N+1)] .dL[1] .tfTarget[1]
%           .xf[5x1] .thrustN .m0kg .ispS .maxIter .tag[char] .kind[char]
% REFERENCES:
%   [1] run_transfer_mee.m:255-258 (res-struct layout); psr_mee_refine.m:294-305
%       (PSR out layout); process/DESIGN_sosc.md sec 4.2.
S = load(matPath);
[~, base] = fileparts(matPath);
geoXf = [1;0;0;0;0];                             % GEO default when fp.xf absent
if isfield(S,'res')                              % MEE_M2 row
    r  = S.res;  fu = r.fuel;  fp = r.fp;
    saved = struct('sigma', r.sigma(:), 'X', fu.X, 'U', fu.U, 'dL', fu.dL, ...
        'tfTarget', r.tf, 'xf', optdef(fp,'xf',geoXf), 'thrustN', fp.thrustN, ...
        'm0kg', fp.m0kg, 'ispS', fp.ispS, 'maxIter', optdef(fp,'maxIter',1500), ...
        'tag', base, 'kind', 'MEE_M2');
elseif isfield(S,'out')                          % PSR-refined row
    o  = S.out;  fu = o.finalOut;  fp = S.fpFinal;
    saved = struct('sigma', o.finalSigma(:), 'X', fu.X, 'U', fu.U, 'dL', fu.dL, ...
        'tfTarget', fp.tf, 'xf', optdef(fp,'xf',geoXf), 'thrustN', fp.thrustN, ...
        'm0kg', fp.m0kg, 'ispS', fp.ispS, ...
        'maxIter', optdef(fp,'maxIter',1500), 'tag', base, 'kind', 'PSR');
else
    error('sosc_load_row:unknownShape', ...
        '%s has neither a res nor an out variable', matPath);
end
end
