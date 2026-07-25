function setup_verify_common()
% SETUP_VERIFY_COMMON  Put orbit_transfer/verify_common on the MATLAB path.
% Self-contained: no campaign paths, no CasADi (callers add CasADi themselves).
% OUTPUTS: none (path side effect)
addpath(fileparts(mfilename('fullpath')));
end
