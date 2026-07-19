function [rhat, that, nhat] = rtn_frame(r, v)
% RTN_FRAME  Local orbital RTN triad (radial, transverse, normal) at inertial
% state (r,v). Columns of [rhat that nhat] rotate RTN->inertial.
% INPUTS:  r [3x1] inertial position; v [3x1] inertial velocity
% OUTPUTS: rhat/that/nhat [3x1 each] unit radial/transverse/normal (inertial)
rhat = r / norm(r);
hvec = cross(r, v);
nhat = hvec / norm(hvec);
that = cross(nhat, rhat);
end
