% name -- tfmin
%
%   Usage
%     tfmin(Tmax)
%
%   Description
%     Estimation of the minmum transfer time for given thrust Tmax
%

% Written on Wed Jan 15 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

function tmin = tfmin(Tmax)

tmin = 850/Tmax; 
