% name -- evolution
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

form1r = ''; for i=1:n form1r = [form1r ' %f']; end;
form1w = ''; for i=1:n form1w = [form1w ' %24.16e']; end;
form2 = ''; for i=1:n form2 = [form2 ' %d']; end;
form3r = ''; for i=1:2*n form3r = [form3r ' %f']; end;
form3w = ''; for i=1:2*n form3w = [form3w ' %24.16e']; end;
form4r = ''; for i=1:lpar form4r = [form4r ' %f']; end;
form4w = ''; for i=1:lpar form4w = [form4w ' %24.16e']; end;
form5 = ''; for i=1:lipar form5 = [form5 ' %d']; end;
form6 = ''; for i=1:n+1 form6 = [form6 ' %f']; end;

% modification 2: 29th july 2007
% author: J. Gergaud
% add 1 because xpoint_{n+2}=|u|
form7 = ''; for i=1:2*n+1+m+2 form7 = [form7 ' %f']; end;
% end of the modification

form8 = ''; for i=1:2*n+1+m form8 = [form8 ' %f']; end;

fid = fopen(name);
s = fscanf(fid,' # %s',[1 1]);
if (~strcmp(s,'path'))
  disp('mfmax: bad file format'); else
  free0 = fscanf(fid,[' # FREE0 = ' form2],[1 n]);
  y0 = fscanf(fid,[' # Y0 = ' form3r],[1 2*n]);
  t0 = fscanf(fid,' # T0 = %f',[1 1]);
  tmin = fscanf(fid,' # TF = %f',[1 1]);
  par = fscanf(fid,[' # PAR = ' form4r],[1 lpar]);
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

%Draw evolution of criterium ?
rep = input('Do you want to draw criterium evolution (yes) ? ','s');
if isempty(rep)|strcmp(rep,'yes')|strcmp(rep,'y')
    
    Nl = length(lambda);

    % save file
    if exist('in.dat','file') ~unix('\cp in.dat inaux.dat'); end;
    if exist('out.dat','file') ~unix('\cp out.dat outaux.dat'); end;
    
    % create initial data file
    infile = 'in.dat';
    ipar(1) = 5; %Solver must only integrate BVP with criterium
    ipar(2) 
    for i=1:Nl
    
    fid = fopen(infile,'w');
    
    fprintf(fid,[form1w '\n'],z(:,i));
    fprintf(fid,[form2 '\n'],free0);
    fprintf(fid,[form3w '\n'],y0);
    fprintf(fid,' %f\n',t0);
    fprintf(fid,' %f\n',tmin);
    fprintf(fid,' %f %f\n',[1. lambda(i)]);
    fprintf(fid,[form4w '\n'],par);
    fprintf(fid,[form5 '\n'],ipar);
    fprintf(fid,' %d \n',1); % onnly one integration step
    % Information for zero path tracking : not needed here
    fprintf(fid,' 1e-4\n');
    fprintf(fid,' 0\n 0. 0.\n');
    fclose(fid);
    
    % Call solver
    ~unix('../bin/path');
    
    fid = fopen('outplus.dat');
	s = fscanf(fid,' # %s',[1 1]);
	if (~strcmp(s,'trajectory'))
 	disp('mfmax: bad file format'); else
	fscanf(fid,[' # FREE0 = ' form2],[1 n]);
        fscanf(fid,[' # Y0 = ' form3r],[1 2*n]);
        fscanf(fid,' # T0 = %f',[1 1]);
        fscanf(fid,' # TF = %f',[1 1]);
        fscanf(fid,[' # LAMBDA = %f %f'],[1 2]);
        fscanf(fid,[' # PAR = ' form4r],[1 lpar]);
        fscanf(fid,[' # IPAR = ' form5],[1 lipar]);
        fscanf(fid,' # NIT = %d',[1 1]);
        fscanf(fid,[' # Z = ' form1r],[1 n]);
        fscanf(fid,[' # S = ' form1r],[1 n]);
        fscanf(fid,' # |S| = %f',[1 1]);
        s = fscanf(fid,' # %s',[1 1]);
        if (~strcmp(s,'Data'))
          disp('mfmax: bad file format 2'); else
          clear data;

% modification 2: 29th july 2007
% author: J. Gergaud
% add one because xpoint_{n+2}=|u|		  
          data = fscanf(fid,form7,[2*n+1+m+2 inf]);
% end of the modification	  
	  Jl(i) = data(1+2*n+m+1,2);
	  

% modification 2: 29th july 2007
% author: J. Gergaud	  
	  J1(i) = data(1+2*n+m+2,2);
% end of the modification
	  
	end;
	end;
        fclose(fid);

    end; %for
    
    figure;
    plot(lambda,Jl,'b');
    hold on;
    plot(lambda,J1,'r');
    legend('J_{\lambda}','J_1');
    
    % reload save files if necessary
    if exist('inaux.dat','file') ~unix('\mv inaux.dat in.dat'); end;
    if exist('outaux.dat','file') ~unix('\mv outaux.dat out.dat'); end;
    
end;

%Draw evolution of control for somme lambda ?
rep = input('Do you want to draw evolution of control (yes) ? ','s');
if isempty(rep)|strcmp(rep,'yes')|strcmp(rep,'y')
   
   lambda
   cpt = 1;
   pos = input('First \lambda (1)');
   if isempty(pos) pos=1; end;
   position(cpt) = pos
   cpt = cpt+1;
   
   while (pos >= 0) 
      pos = input('Give next position of \lambda (-1 to stop) ');
      if isempty(pos) pos=cpt; end;
      if (pos > 0)
	 position(cpt) = pos;
	 cpt = cpt+1;
      end;
   end; % while
   
   Npos = cpt-1;
   u = zeros(Npos,m,nit+1);
   norm_u = zeros(Npos,nit+1);
   L = zeros(Npos,nit+1);
   
   % create initial data file
   infile = 'in.dat';
   ipar(1) = 3; %Solver must only integrate BVP

   for i = 1:Npos
     fid = fopen(infile,'w');
    
     fprintf(fid,[form1w '\n'],z(:,position(i)));
     fprintf(fid,[form2 '\n'],free0);
     fprintf(fid,[form3w '\n'],y0);
     fprintf(fid,' %f\n',t0);
     fprintf(fid,' %f\n',tmin);
     fprintf(fid,' %f %f\n',[1. lambda(position(i))]);
     fprintf(fid,[form4w '\n'],par);
     fprintf(fid,[form5 '\n'],ipar);
     fprintf(fid,' %d \n',nit);
     % Information for zero path tracking : not needed here
     fprintf(fid,' 1e-4\n');
     fprintf(fid,' 0\n 0. 0.\n');
     fclose(fid);
    
     % Call solver
     ~unix('../bin/path');
     
     % Read data
     fid = fopen('out.dat');
     s = fscanf(fid,' # %s',[1 1]);
     if (~strcmp(s,'trajectory'))
       disp('mfmax: bad file format'); else
       fscanf(fid,[' # FREE0 = ' form2],[1 n]);
       fscanf(fid,[' # Y0 = ' form3r],[1 2*n]);
       fscanf(fid,' # T0 = %f',[1 1]);
       fscanf(fid,' # TF = %f',[1 1]);
       fscanf(fid,[' # LAMBDA = %f %f'],[1 2]);
       fscanf(fid,[' # PAR = ' form4r],[1 lpar]);
       fscanf(fid,[' # IPAR = ' form5],[1 lipar]);
       fscanf(fid,' # NIT = %d',[1 1]);
       fscanf(fid,[' # Z = ' form1r],[1 n]);
       fscanf(fid,[' # S = ' form1r],[1 n]);
       fscanf(fid,' # |S| = %f',[1 1]);
       s = fscanf(fid,' # %s',[1 1]);
       if (~strcmp(s,'Data'))
         disp('mfmax: bad file format 2'); else
         clear data;
         data = fscanf(fid,form8,[2*n+1+m inf]);
	  
         u(i,:,:) = data(2*n+2:2*n+1+m,:);
         norm_u(i,:) = sqrt(u(i,1,:).^2 + u(i,2,:).^2 + u(i,3,:).^2);
         L(i,:) = data(n,:);
       end;
     end;
     fclose(fid);

   end; % for
   
   t = data(1,:);
   
   % color reserve
% modofication: 28th july 2007
% author: J. Gergaud
% modification of the order of the colors
   col(1) = 'b';
   col(2) = 'c';
   col(3) = 'k';
   col(4) = 'g';
   col(5) = 'r';
% end of the modification

   figure;
   hold on;
   for i = 1:Npos
% modification 28th july 2007
% author: J. Gergaud   
      j = mod(i-1,5)+1;
%     j = mod(i,5)+1;
      plot(t,norm_u(i,:),col(j));
%	  plot(t,norm_u(i,:),col(mod(i,6)));
% end of the modification
   end;
   
end; % evolution control

end;
end;
