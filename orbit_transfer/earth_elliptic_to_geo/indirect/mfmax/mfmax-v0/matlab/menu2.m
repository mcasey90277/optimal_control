% name -- menu2
%
%  Usage
%    menu2
%
%  Description
%    call solver
%
%  See also
%    menu1, menu3, menu4
%

% Written on Thu Nov 27 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

r = 0; while ~r
  infile = input('Input file name (in.dat) ? ','s');
  if isempty(infile) infile = 'in.dat'; end;
  r = exist(infile,'file');
  if ~r
    disp(['mfmax: ' infile ' not found']);
  else
    if ~strcmp(infile,'in.dat')
      r = ~unix(['\cp ' infile ' in.dat']);
    end;
  end;
end;

~unix(['../bin/path']);

% Read infile to know which are output files

n = 7;
lpar = 14;
lipar = 2;

form1 = ''; for i = 1:n form1 = [ form1 ' %f' ]; end;
form2 = ''; for i = 1:n form2 = [ form2 ' %d' ]; end;
form3 = ''; for i = 1:(2*n) form3 = [ form3 ' %f' ]; end;
form4 = ''; for i = 1:lpar form4 = [ form4 ' %f' ]; end;
form5 = ''; for i = 1:lipar form5 = [ form5 ' %d']; end;

fid = fopen('in.dat');
zi = fscanf(fid,form1,[1 n]);
free0 = fscanf(fid,form2,[1 n]);
y0 = fscanf(fid,form3,[1 2*n]);
L0 = fscanf(fid,' %f', [1 1]);
Lmin = fscanf(fid,' %f',[1 1]);
lambda = fscanf(fid,' %f %f',[1 2]);
par = fscanf(fid,form4,[1 lpar]);
ipar = fscanf(fid,form5,[1 lipar]);
nit = fscanf(fid,' %d',[1 1]);
jacstep = fscanf(fid,' %f',[1 1]);
trace = fscanf(fid,' %d',[1 1]);
fclose(fid);

if (ipar(1) <= 2)
    outext = input('Solution and trace file extension ()? ','s');
    outfile = ['out' outext '.dat'];
    tracefile = ['path' outext '.trc'];
    file = ['fort.' num2str(trace)];
    if ~strcmp(outfile,'out.dat')
        ~unix(['\mv out.dat ' outfile]);
        ~unix(['\mv ' file ' ' tracefile]);
    else
        ~unix(['\mv ' file ' ' tracefile]);
    end;
else
    outext = input('Solution file extension ()? ','s');
    outfile = ['out' outext '.dat'];
    if ~strcmp(outfile,'out.dat')
        ~unix(['\mv out.dat ' outfile]);
    end;
end;
