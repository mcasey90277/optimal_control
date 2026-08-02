% name -- drawpath
%
%   Usage
%     drawpath(name)
%
%   Description
%     Draw zero path
%

% Written on Thu Nov 27 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

function drawpath(name)

format long e;
n = 7;
m = 3;
lpar = 14;
lipar = 2;
   
form1 = ''; for i=1:n form1 = [form1 ' %f']; end;
form1w = ''; for i=1:n form1w = [form1w ' %24.16e']; end;
form2 = ''; for i=1:n form2 = [form2 ' %d']; end;
form3 = ''; for i=1:2*n form3 = [form3 ' %f']; end;
form3w = ''; for i=1:2*n form3w = [form3w ' %24.16e']; end;
form4 = ''; for i=1:lpar form4 = [form4 ' %f']; end;
form4w = ''; for i=1:lpar form4w = [form4w ' %24.16e']; end;
form5 = ''; for i=1:lipar form5 = [form5 ' %d']; end;
form6 = ''; for i=1:n+1 form6 = [form6 ' %f']; end;

fid = fopen(name);
s = fscanf(fid,' # %s',[1 1]);
if (~strcmp(s,'path'))
  disp('mfmax: bad file format'); else
  free0 = fscanf(fid,[' # FREE0 = ' form2],[1 n]);
  y0 = fscanf(fid,[' # Y0 = ' form3],[1 2*n]);
  t0 = fscanf(fid,' # T0 = %f',[1 1]);
  tmin = fscanf(fid,' # TF = %f',[1 1]);
  par = fscanf(fid,[' # PAR = ' form4],[1 lpar]);
  ipar = fscanf(fid,[' # IPAR = ' form5],[1 lipar]);
  nit = fscanf(fid,[' # NIT = %d'],[1 1]);
  s = fscanf(fid,' # %s',[1 1]);
  if (~strcmp(s,'Data'))
    disp('mfmax: bad file format 2'); else
  
    data = fscanf(fid,form6,[n+1 inf]);
fclose(fid);

Tmax = par(5);
Beta = par(6);
mu0 = par(7);
tf = (tmin-t0)*par(1)+t0;

clear z;
z = data(2:n+1,:);
lambda =  data(1,:);

% Draw zero path
rep = input('Draw zero path (yes)? ','s');
if isempty(rep)|strcmp(rep,'y')|strcmp(rep,'yes');
figure;
subplot(7,1,1);
hold on;
title(['Zero path for T_{max} = ' num2str(Tmax) 'N, t_f = ' num2str(tf) 'rad.']);
plot(z(1,:),lambda);
ax = axis;
axis([ax(1) ax(2) 0. 1.1]);

subplot(7,1,2);
plot(z(2,:),lambda);
ax = axis;
axis([ax(1) ax(2) 0. 1.1]);

subplot(7,1,3);
plot(z(3,:),lambda);
ax = axis;
axis([ax(1) ax(2) 0. 1.1]);

subplot(7,1,4);
plot(z(4,:),lambda);
ax = axis;
axis([ax(1) ax(2) 0. 1.1]);

subplot(7,1,5);
plot(z(5,:),lambda);
ax = axis;
axis([ax(1) ax(2) 0. 1.1]);

subplot(7,1,6);
plot(z(6,:),lambda);
ax = axis;
axis([ax(1) ax(2) 0. 1.1]);

subplot(7,1,7);
plot(z(7,:),lambda);
ax = axis;
axis([ax(1) ax(2) 0. 1.1]);

end;


%Draw solution for a specific lambda ?
rep = input('Do you want to draw figures for a specific \lambda (no) ? ','s');
if strcmp(rep,'yes')|strcmp(rep,'y')
    
    lambda
    Nl = length(lambda);
    pos = -1;
    while (pos < 1)|(pos > Nl)
        pos = input('Give position of the \lambda you wish (1) ');
        if isempty(pos) pos=1; end;
        if (pos < 1)|(pos > Nl)
            disp('mfmax: there is no \lambda at this position');
        end;
    end;
    
    % save file
    if exist('in.dat','file') ~unix('\cp in.dat inaux.dat'); end;
    if exist('out.dat','file') ~unix('\cp out.dat outaux.dat'); end;
    
    % create initial data file
    infile = 'in.dat';
    ipar(1) = 3; %Solver must only integrate BVP
    fid = fopen(infile,'w');
    
    fprintf(fid,[form1w '\n'],z(:,pos));
    fprintf(fid,[form2 '\n'],free0);
    fprintf(fid,[form3w '\n'],y0);
    fprintf(fid,' %f\n',t0);
    fprintf(fid,' %f\n',tmin);
    fprintf(fid,' %f %f\n',[1. lambda(pos)]);
    fprintf(fid,[form4w '\n'],par);
    fprintf(fid,[form5 '\n'],ipar);
    fprintf(fid,' %d \n',nit);
    % Information for zero path tracking : not needed here
    fprintf(fid,' 1e-4\n');
    fprintf(fid,' 0\n 0. 0.\n');
    fclose(fid);
    
    % Call solver
    ~unix('../bin/path');
    
    drawres('out.dat');
    
    % reload save files if necessary
    if exist('inaux.dat','file') ~unix('\mv inaux.dat in.dat'); end;
    if exist('outaux.dat','file') ~unix('\mv outaux.dat out.dat'); end;
    
end;

end;
end;
