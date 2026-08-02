% name -- menu3
%
%  Usage
%    menu3
%
%  Description
%    draw some graph
%
%  See also
%    menu1, menu2
%

% Written on Mon Jan 20 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

tr = 1;
r = 0; while ~r
  extension = input('Solution and trace file extension ? ','s');
  outfile = ['out' extension '.dat'];
  tracefile = ['path' extension '.trc'];
  
  r = exist(outfile,'file');
  if ~r
     disp(['mfmax: ' outfile ' not found']);
  end;
  r = exist(tracefile,'file');
  if ~r
     disp(['mfmax: ' tracefile ' not found']);
     rep = input('Do you want to continu anyway (no)? ','s');
     if strcmp(rep,'yes')|strcmp(rep,'y')
        r = true; tr = 0;
     end;
   end;
end;
   
rep = input('Extract graphs from solution file (yes)? ','s');
if isempty(rep)|strcmp(rep,'y')|strcmp(rep,'yes')
   drawres(outfile);
end;
if (tr ~= 0)
   rep = input('Extract graphs from trace file (yes) ? ','s');
   if isempty(rep)|strcmp(rep,'y')|strcmp(rep,'yes')
      drawpath(tracefile);
   end;
end;
   
