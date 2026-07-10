function Z = sms_pack(lam0, yInterior)
% SMS_PACK  Stack Sundman-domain multiple-shooting unknowns into one vector.
%
% INPUTS:
%   lam0      - initial costates [8x1] ([lamR; lamV; lamM; lamT])
%   yInterior - augmented 16-dim states at interior arc joints [16x(M-1)]
%
% OUTPUTS:
%   Z - unknown vector [(16M-8)x1]: [lam0; yInterior(:)]

Z = [lam0(:); yInterior(:)];
end
