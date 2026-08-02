% name -- mfmax
%
%  Usage
%    mfmax
%
%  Description
%    Matlab interface for the maximum final mass transfer
%    problem (3D model, fixed final longitude).
%
%  See also
%    menu1, menu2, menu3, menu4
%

% Written on Tue Apr 10 16:50:00 MET DST 2001
% by Jean-Baptiste Caillau - ENSEEIHT-IRIT, UMR CNRS 5505
% Modified on Fri Jan 17 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

title = 'mfmax - Maximum final Mass transfer';

entries = ['Create Initial data        '; ...
           'Call solver                '; ...
	   'Make graphs                '; ...
	   'Demo                       '; ...
           'Help                       '];

lentries = size(entries,1);

ch = -1; while (ch ~= 0)
% Display main menu
fprintf(1,'\n');
fprintf(1,'\t\t%s\n',title);
fprintf(1,'\n');
for i = 1:lentries
  fprintf(1,'\t\t%d. %s\n',i,entries(i,:));
end;
fprintf(1,'\n');
fprintf(1,'\t\t0. Quit\n');
fprintf(1,'\n');

% Choice
ch = -1; while (ch == -1)
  fprintf(1,'Choice ? '); ch = input('');
  if isempty(ch) ch = -1; end;
  if (ch < 0) | (ch > lentries)
    ch = -1;
  elseif (ch == 0)
    close all;
  else
    fprintf(1,'\n%d. %s\n\n',ch,entries(ch,:));
    eval(['menu' num2str(ch)]);
  end;
end;

end;

