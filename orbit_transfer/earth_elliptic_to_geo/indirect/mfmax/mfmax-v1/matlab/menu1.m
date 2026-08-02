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

% Written on Thu Nov 27 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505


format long e;
n = 8;

% Maximum thrust
Tmax_d = 10.0;
fprintf(1,'**** Enter the maximum thrust ****\n');
fprintf(1,'Tmax in Newtons (%f) ? ',Tmax_d); Tmax = input(''); 
if isempty(Tmax) Tmax = Tmax_d; end;

fprintf(1,'Do you want default values for problem definition (yes) ? ');
rep = input('');

% Default values
Beta = 1.42e-2;
mu0  = 5165.8620912;
P0_d = 11.625;
ex0_d = 0.75;
ey0_d = 0.0;
hx0_d = 0.0612;
hy0_d = 0.0;
m0_d = 1500; 
Pf_d = 42.165;
exf_d = 0.0;
eyf_d = 0.0;
hxf_d = 0.0;
hyf_d = 0.0;
nit_d = min(floor(5000*clf/Tmax),8000);

if isempty(rep)|strcmp(rep,'y')|strcmp(rep,'yes')
    P0  = P0_d;
    ex0 = ex0_d;
    ey0 = ey0_d;
    hx0 = hx0_d;
    hy0 = hy0_d;
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
    fprintf(1,'P0 in Megameters(%f) ? ',P0_d); P0 = input('');
    if isempty(P0) P0 = P0_d; end;
    fprintf(1,'ex0 (%f) ? ',ex0_d); ex0 = input(''); 
    if isempty(ex0) ex0 = ex0_d; end;
    fprintf(1,'ey0 (%f) ? ',ey0_d); ey0 = input(''); 
    if isempty(ey0) ey0 = ey0_d; end;
    fprintf(1,'hx0 (%f) ? ',hx0_d); hx0 = input(''); 
    if isempty(hx0) hx0 = hx0_d; end;
    fprintf(1,'hy0 (%f) ? ',hy0_d); hy0 = input(''); 
    if isempty(hy0) hy0 = hy0_d; end;
    fprintf(1,'m0 in kilogrammes (%f) ? ',m0_d); m0 = input(''); 
    if isempty(m0) m0 = m0_d; end;
    fprintf(1,'Pf in Megameters (%f) ? ',Pf_d); Pf = input(''); 
    if isempty(Pf) Pf = Pf_d; end;
    fprintf(1,'exf (%f) ? ',exf_d); exf = input(''); 
    if isempty(exf) exf = exf_d; end;
    fprintf(1,'eyf (%f) ? ',eyf_d); eyf = input(''); 
    if isempty(eyf) eyf = eyf_d; end;
    fprintf(1,'hxf (%f) ? ',hxf_d); hxf = input(''); 
    if isempty(hxf) hxf = hxf_d; end;
    fprintf(1,'hyf (%f) ? ',hyf_d); hyf = input(''); 
    if isempty(hyf) hyf = hyf_d; end;
    
    fprintf(1,'nit (%d) ? ',nit_d); nit = input('');
    if isempty(nit) nit = nit_d; end;

end;

fprintf(1,'\n');
fprintf(1,'**** Summary of main physical values in standard units ****\n');
fprintf(1,'Tmax\t= %f Newtons\n',Tmax);
fprintf(1,'P0\t= %f Mmeters\n',P0);
fprintf(1,'ex0\t= %f\n',ex0);
fprintf(1,'ey0\t= %f\n',ey0);
fprintf(1,'hx0\t= %f\n',hx0);
fprintf(1,'hy0\t= %f\n',hy0);
fprintf(1,'m0\t= %f kilogrammes\n',m0);
fprintf(1,'Pf\t= %f Mmeters\n',Pf);
fprintf(1,'exf\t= %f\n',exf);
fprintf(1,'eyf\t= %f\n',eyf);
fprintf(1,'hxf\t= %f\n',hxf);
fprintf(1,'hyf\t= %f\n',hyf);

% Final Longitude
fprintf(1,'Initial longitude (%f) ?\n',pi);
L0 = input('');
if isempty(L0)  L0 = pi; end;
fprintf(1,'\n**** Enter minimum transfert longitude ****\n');
lmin_d = Lfmin(Tmax,L0);
fprintf(1,'lfmin in radian (%f rad.) ? ',lmin_d);
Lmin = input('');
if isempty(Lmin) Lmin = lmin_d; end;
clf_d = 2.;
fprintf(1,'Which multiplier for Lf (%f) ? ', clf_d);
clf = input('');
if isempty(clf) clf = clf_d; end;

% Step for initial condition discrete continuation
rep = input('Do you want default values for discrete continuation (yes)?','s');
if isempty(rep)|strcmp(rep,'yes')|strcmp(rep,'y')
   stepci = 0.1;
   minstepci = 0.01;
else
   stepci_d = 0.1;
   minstepci_d = 0.01;
   fprintf(1,'Initial and maximum (same) step (%f)?\n',stepci_d);
   stepci = input('');
   if ispempty(stepci) stepci = stepci_d; end;
   fprintf(1,'Minimum step (%f)?\n',minstepci_d);
   minstepci = input('');
   if isempty(minstepci) minstepci = minstepci_d; end;
end;

zi = zeros(1,n);
zi(1) = 10.; % zi(1) cannot be <= 0
free0 = [8 9 10 11 12 13 14 15];
y0 = [P0 ex0 ey0 hx0 hy0 0. m0 0. 0. 0. 0. 0. 0. 0. 0. 0.];
lambda = [0. 0.];
par = [clf Tmax Beta mu0 Pf exf eyf hxf hyf stepci minstepci];
lpar = length(par);
lipar = 1;
ipar(1) = 1;

% Option file
rep = input('\nDo you want default parameters for zero curve tracking (yes)? ','s');
if isempty(rep)|strcmp(rep,'yes')|strcmp(rep,'y')
    jac_step = 1e-5;
    trace = 2;
    sspar = [0. 0.];
else
    jac_step_d = 1e-5;
    trace = 2;
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
form2 = ''; for i = 1:(2*n) form2 = [ form2 '%15.15e ' ]; end;
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
fprintf(fid,'%15.15e\n',L0);
fprintf(fid,'%15.15e\n',Lmin);
fprintf(fid,' %15.15e   %15.15e\n',lambda);
fprintf(fid,[form3 '\n'],par);
fprintf(fid,'  %d \n',ipar);
fprintf(fid,'%d\n',nit);
fprintf(fid,'  %15.15e\n',jac_step);
fprintf(fid,'  %d\n',trace);
fprintf(fid,'  %15.15e   %15.15e\n',sspar);

fclose(fid);
