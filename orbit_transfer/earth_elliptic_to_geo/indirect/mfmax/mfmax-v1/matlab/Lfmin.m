% name -- Lfmin
%
%   Usage
%     Lfmin(L0,Tmax)
%
%   Description
%     Estimation of the minmum transfer longitude for given thrust Tmax
%

% Written on Wed Jan 15 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

function Lmin = Lfmin(Tmax,L0)

ref = 267.537780225 ;

Lmin = ref/Tmax + L0 ;
