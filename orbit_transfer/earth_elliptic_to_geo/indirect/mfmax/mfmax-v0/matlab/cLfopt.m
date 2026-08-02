% name -- cLfopt
%
%   Usage
%     cLfopt(ctf,Tmax)
%
%   Description
%     Estimation of optimal Lf multiplier for givrn ctf & Tmax
%

% Written on Wed Jan 15 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

function cLopt = cLfopt(Tmax,ctf)

cLopt = 1.12*ctf + 0.09;
 
