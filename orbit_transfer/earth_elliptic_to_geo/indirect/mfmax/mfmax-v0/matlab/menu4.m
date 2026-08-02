% name -- menu4
%
%  Usage
%    menu4
%
%  Description
%    Call the solver for the demo input files and draw figures
%
%  See also
%    menu1, menu2, menu3, menu5
%

% Written on Thu Nov 27
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

infile = 'IN/in10N.dat'

if ~exist(infile)
  disp(['mfmax: ' infile ' not found']);
else
  r = ~unix(['\cp ' infile ' in.dat']);
  if ~r
    disp('mfmax: cannot copy input file');
  else
    r = ~unix('../bin/path');
    if ~r disp('mfmax: cannot execute mfmax'); end;
  end;
end;

if r
  unix('\mv fort.2 path.trc');
  drawpath('path.trc');
  drawres('out.dat');
end;
