% name -- menu1
%
%  Usage
%    menu1
%
%  Description
%    create initial data
%
%  See also
%    menu2, menu3, menu4, menu5
%

% Written on Fri Jan 17 2003
% Last modified on Wed Jan 7 2004
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505


format long e;
n = 7;

% Maximum thrust
Tmax_d = 10.0;
fprintf(1,'**** Enter the maximum thrust ****\n');
fprintf(1,'Tmax in Newtons (%15.15e) ? ',Tmax_d); Tmax = input(''); 
if isempty(Tmax) Tmax = Tmax_d; end;
% Move Tmax into good unit
Tmax = Tmax;

fprintf(1,'Do you want default values for problem definition (enter for yes) ? ');
rep = input('');

% Default values
Beta = 1.42e-2;
mu0  = 5165.8620912;
P0_d = 11.625;
ex0_d = 0.75;
ey0_d = 0.0;
hx0_d = 0.0612;
hy0_d = 0.0;
L0_d = 3.14159265358979;
m0_d = 1500; 
Pf_d = 42.165;
exf_d = 0.0;
eyf_d = 0.0;
hxf_d = 0.0;
hyf_d = 0.0;
nit_d = 8000;

if isempty(rep)
    P0  = P0_d;
    ex0 = ex0_d;
    ey0 = ey0_d;
    hx0 = hx0_d;
    hy0 = hy0_d;
    L0  = L0_d;
    m0  = m0_d;
    
    Pf  = Pf_d;
    exf = exf_d;
    eyf = eyf_d;
    hxf = hxf_d;
    hyf = hyf_d;
    nit = nit_d;
        
else
    
    
    fprintf(1,'\n');
    fprintf(1,'**** Enter values ****\n');
    fprintf(1,'P0 in Megameters(%15.15e) ? ',P0_d); P0 = input('');
    if isempty(P0) P0 = P0_d; end;
    fprintf(1,'ex0 (%15.15e) ? ',ex0_d); ex0 = input(''); 
    if isempty(ex0) ex0 = ex0_d; end;
    fprintf(1,'ey0 (%15.15e) ? ',ey0_d); ey0 = input(''); 
    if isempty(ey0) ey0 = ey0_d; end;
    fprintf(1,'hx0 (%15.15e) ? ',hx0_d); hx0 = input(''); 
    if isempty(hx0) hx0 = hx0_d; end;
    fprintf(1,'hy0 (%15.15e) ? ',hy0_d); hy0 = input(''); 
    if isempty(hy0) hy0 = hy0_d; end;
    fprintf(1,'L0 in radians (%15.15e) ? ',L0_d); L0 = input(''); 
    if isempty(L0) L0 = L0_d; end;
    fprintf(1,'m0 in kilogrammes (%15.15e) ? ',m0_d); m0 = input(''); 
    if isempty(m0) m0 = m0_d; end;
    fprintf(1,'Pf in Megameters (%15.15e) ? ',Pf_d); Pf = input(''); 
    if isempty(Pf) Pf = Pf_d; end;
    fprintf(1,'exf (%15.15e) ? ',exf_d); exf = input(''); 
    if isempty(exf) exf = exf_d; end;
    fprintf(1,'eyf (%15.15e) ? ',eyf_d); eyf = input(''); 
    if isempty(eyf) eyf = eyf_d; end;
    fprintf(1,'hxf (%15.15e) ? ',hxf_d); hxf = input(''); 
    if isempty(hxf) hxf = hxf_d; end;
    fprintf(1,'hyf (%15.15e) ? ',hyf_d); hyf = input(''); 
    if isempty(hyf) hyf = hyf_d; end;
    
    fprintf(1,'nit (%15.15e) ? ',nit_d); nit = input('');
    if isempty(nit) nit = nit_d; end;

end;

fprintf(1,'\n');
fprintf(1,'**** Summary of main physical values in standard units ****\n');
fprintf(1,'Tmax\t= %15.15e Newtons\n',Tmax);
fprintf(1,'P0\t= %15.15e Mmeters\n',P0);
fprintf(1,'ex0\t= %15.15e\n',ex0);
fprintf(1,'ey0\t= %15.15e\n',ey0);
fprintf(1,'hx0\t= %15.15e\n',hx0);
fprintf(1,'hy0\t= %15.15e\n',hy0);
fprintf(1,'L0\t= %15.15e radians\n',L0);
fprintf(1,'m0\t= %15.15e kilogrammes\n',m0);
fprintf(1,'Pf\t= %15.15e Mmeters\n',Pf);
fprintf(1,'exf\t= %15.15e\n',exf);
fprintf(1,'eyf\t= %15.15e\n',eyf);
fprintf(1,'hxf\t= %15.15e\n',hxf);
fprintf(1,'hyf\t= %15.15e\n',hyf);

% Final time
tmin_d = tfmin(Tmax);
fprintf(1,'\n**** Enter minimum transfert time ****\n');
fprintf(1,'tfmin in hours (%f h.) ? ',tmin_d);
tmin = input('');
if isempty(tmin) tmin = tmin_d; end;
ctf_d = 1.5;
fprintf(1,'Which multiplier for tf (%f) ? ', ctf_d); ctf = input('');
if isempty(ctf) ctf = ctf_d; end;

% Final longitude
Lmin_d = Lfmin(Tmax,L0);
fprintf(1,'\n**** Enter data for final longitude ****\n');
fprintf(1,'Lfmin in radian is (%f) ? ',Lmin_d); Lmin = input('');
if isempty(Lmin) Lmin = Lmin_d; end;
cLf_d = cLfopt(Tmax,ctf);
rep = input('Solve for a given Lf (yes)? ','s');
if isempty(rep)|strcmp(rep,'yes')|strcmp(rep,'y')
     fprintf(1,'Which multiplier for Lfmin (%f)? ',cLf_d);
     cLf = input('');
     if isempty(cLf) cLf = cLf_d; end;
     Lffix = 1;
     tolLf = 0.01;
else
     fprintf(1,'Approximation of Lf optimal multiplier (%f)? ',cLf_d);
     cLf = input('');
     if isempty(cLf) cLf = cLf_d; end;
     tolLf_d = 0.01;
     fprintf(1,'Error tolerance for Lf optimal multiplier (%f)? ',tolLf_d);         tolLf = input('');
     if isempty(tolLf) tolLf = tolLf_d; end;
     Lffix = 2;
end;

% steps for homotopy on initial condition
if (Tmax >= 1.)
  stepci_d = 0.2;
else
  stepci_d = 0.1;
end;
fprintf(1,'Maximum (and starting) step for initial condition homotopy (%f)? ',stepci_d);
stepci = input('');
if isempty(stepci) stepci = stepci_d; end;
stepcimin_d = 1.e-2;
fprintf(1,'Minimum step for initial condition homotopy (%f)? ',stepcimin_d);
stepcimin = input('');
if isempty(stepcimin) stepcimin = stepcimin_d; end;

% Homotopic criterium
rep = input('\nConvex (c) or power (p) criterium (p)? ','s');
if isempty(rep) crit = 2;
else   crit = 1;
end;

zi = zeros(1,n);
free0 = [8 9 10 11 12 13 14];
y0 = [P0 ex0 ey0 hx0 hy0 L0 m0 0. 0. 0. 0. 0. 0. 0.];
t0 = 0.;
tf = tmin;
lambda = [0. 0.];
par = [ctf tolLf cLf Lmin Tmax Beta mu0 stepci stepcimin Pf exf eyf hxf hyf];
lpar = length(par);
ipar = [Lffix crit];
lipar = 2;

% Option file
rep = input('\nDo you want default parameters for zero curve tracking (yes)? ','s');
if isempty(rep)|strcmp(rep,'yes')|strcmp(rep,'y')
    jac_step = 1e-4;
    trace = 1;
    sspar = [0. 0.];
else
    jac_step_d = 1e-4;
    trace = 1;
    sspar_d = [0. 0.];
    sspar = sspar_d;
    fprintf(1,'Which step for jacobian evaluation (%15.15e) ? ',jac_step_d);jac_step = input('');
    if isempty(jac_step) jac_step = jac_step_d; end;
    fprintf(1,'Minimum step size (small) ? '); sspar(1) = input('');
    if isempty(sspar(1)) sspar(1) = 0.; end;
    fprintf(1,'Maximum step size (1.) ? '); sspar(2) = input('');
    if isempty(sspar(2)) sspar(2) = 0.; end;
        
end;

form1 = ''; for i = 1:n form1 = [ form1 '%15.15e ' ]; end;
form2 = ''; for i = 1:(2*n) form2 = [ form2 '%15.12e ' ]; end;
form3 = ''; for i = 1:lpar form3 = [ form3 '%15.15e ' ]; end;
form4 = ''; for i = 1:n form4 = [ form4 '%d ' ]; end;

fprintf(1,'\n');
fid = -1; while (fid == -1) 
    infile = input('Input file name (in.dat) ? ','s');
    if isempty(infile) infile = 'in.dat'; end;
    fid = fopen(infile,'w'); 
    if (fid == -1) disp('maxmf: unable to open file'); end;
end;

fprintf(fid,[form1 '\n'],zi);
fprintf(fid,[form4 '\n'],free0);
fprintf(fid,[form2 '\n'],y0);
fprintf(fid,'%15.15e\n',t0);
fprintf(fid,'%15.15e\n',tf);
fprintf(fid,' %15.15e   %15.15e\n',lambda);
fprintf(fid,[form3 '\n'],par);
fprintf(fid,'  %d  %d \n',ipar);
fprintf(fid,'%d\n',nit);
fprintf(fid,'  %15.15e\n',jac_step);
fprintf(fid,'  %d\n',trace);
fprintf(fid,'  %15.15e   %15.15e\n',sspar);

fclose(fid);
