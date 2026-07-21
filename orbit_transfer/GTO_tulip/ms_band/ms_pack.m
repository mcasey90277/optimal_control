function Z = ms_pack(lam0, yInterior)
% MS_PACK  Stack multiple-shooting unknowns into a single vector.
%
% INPUTS:
%   lam0      - initial costates [7x1]
%   yInterior - augmented states at interior arc joints [14x(M-1)]
%
% OUTPUTS:
%   Z - unknown vector [(14M-7)x1]: [lam0; yInterior(:)]

Z = [lam0(:); yInterior(:)];
end
