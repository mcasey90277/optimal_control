function tf = same_root(zA, zB)
%% Purpose:
%
%   Are two solution vectors the same root? Both must be finite with
%   non-zero costates; the seven initial costates must agree to 1e-6 of the
%   larger norm (SYMMETRIC), and the flight times to 1e-6 relative.
%
%  ASSUMPTIONS / NOTES:
%
% • A re-polish moves a root by ~1e-9; this is a registry rule, not a proof
%   of identity -- near a fold distinct roots can be closer than any fixed
%   tolerance, and an ill-conditioned polish can move one root further
%   (review 2026-09-19). Normal chart, so there is no scale to quotient.
%
%% Inputs:
%
%  zA, zB                   double [8 x 1]          costates 1:7, t_f 8
%
%% Outputs:
%
%  tf                       logical                 same root
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

tf = false;
if numel(zA) < 8 || numel(zB) < 8 || ~all(isfinite(zA(1:8))) || ~all(isfinite(zB(1:8))), return, end
a = zA(1:7);  b = zB(1:7);
na = sqrt(sum(a(:).^2));  nb = sqrt(sum(b(:).^2));
if na == 0 || nb == 0, return, end
tf = sqrt(sum((a(:) - b(:)).^2)) <= 1e-6*max(na, nb) && abs(zA(8) - zB(8)) <= 1e-6*max(abs(zA(8)), abs(zB(8)));
end
