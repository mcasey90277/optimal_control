% name -- menu5
%
%  Usage
%    menu5
%
%  Description
%    open the refence guide with acrobat
%
%  See also
%    menu1, menu2, menu3, menu4
%

acrobat = 'acrobat';
dir = '../doc';
name = 'mfmax.pdf';
r = 0; while ~r
  r = ~unix(sprintf('%s %s/%s',acrobat,dir,name));
  if ~r
    fprintf(1,'mfmax: cannot execute %s\n',acrobat);
    acrobat = input('Name of your PDF reader (cancel) ? ','s');
    if isempty(acrobat) r = 1; end;
  end;
end;

% Written on Tue Apr 10 16:50:00 MET DST 2001
% by Jean-Baptiste Caillau - ENSEEIHT-IRIT, UMR CNRS 5505
% Very slightly modified on Thu Nov 27 14:53:36 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505
