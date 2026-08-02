% name -- drawres
%
%   Usage
%     drawres(name)
%
%   Description
%     Draw solution
%

% Written on Thu Nov 27 2003
% by Thomas Haberkorn - ENSEEIHT-IRIT, UMR CNRS 5505

function drawres(name)
     
close all;

n = 7;
m = 3;
lpar = 14;
lipar = 2;
   
form1 = ''; for i=1:n form1 = [form1 ' %f']; end;
form2 = ''; for i=1:n form2 = [form2 ' %d']; end;
form3 = ''; for i=1:2*n form3 = [form3 ' %f']; end;
form4 = ''; for i=1:lpar form4 = [form4 ' %f']; end;
form5 = ''; for i=1:lipar form5 = [form5 ' %d']; end;
form6 = ''; for i=1:2*n+1+m form6 = [form6 ' %f']; end;

fid = fopen(name);
s = fscanf(fid,' # %s',[1 1]);
if (~strcmp(s,'trajectory'))
  disp('mfmax: bad file format'); else
  free0 = fscanf(fid,[' # FREE0 = ' form2],[1 n]);
  y0 = fscanf(fid,[' # Y0 = ' form3],[1 2*n]);
  t0 = fscanf(fid,' # T0 = %f',[1 1]);
  tmin = fscanf(fid,' # TF = %f',[1 1]);
  lambda = fscanf(fid,[' # LAMBDA = %f %f'],[1 2]);
  par = fscanf(fid,[' # PAR = ' form4],[1 lpar]);
  ipar = fscanf(fid,[' # IPAR = ' form5],[1 lipar]);
  nit = fscanf(fid,' # NIT = %d',[1 1]);
  z = fscanf(fid,[' # Z = ' form1],[1 n]);
  S = fscanf(fid,[' # S = ' form1],[1 n]);
  norm_S = fscanf(fid,' # |S| = %f',[1 1]);
  s = fscanf(fid,' # %s',[1 1]);
  if (~strcmp(s,'Data'))
    disp('mfmax: bad file format 2'); else
    clear data;
    data = fscanf(fid,form6,[2*n+1+m inf]);
fclose(fid);

Tmax = par(5);
Beta = par(6);
mu0 = par(7);

t = data(1,:); 
N = length(t);
tf = data(1,N);
P  = data(2,:);
ex = data(3,:);
ey = data(4,:);
hx = data(5,:);
hy = data(6,:);
L  = data(7,:);
M  = data(8,:);
p  = data(9:2*n+1,:);
u  = data(2*n+2:2*n+1+m,:);

norm_u = sqrt(u(1,:).^2 + u(2,:).^2 + u(3,:).^2) ;

Lm  = (mod(L,2*pi)/pi)-1;

fig1 = 1;
fig2 = 2;
fig3 = 3;

% Draw graph with only control
figure(fig1);
subplot(4,1,1);
title(['Control vs time, T_{max} = ' num2str(Tmax) ' N, ctf = ' num2str(par(1)) ...
       ', m_f = ' num2str(M(N)) ' kg.']);
hold on;
plot(t,u(1,:));
xlabel('time');
ylabel('radial thrust');
ax = axis;
axis([t(1) t(N) -1.1 1.1]);

subplot(4,1,2);
hold on;
plot(t,u(2,:));
xlabel('time');
ylabel('ortho-radial thrust');
ax = axis;
axis([t(1) t(N) -1.1 1.1]);

subplot(4,1,3);
hold on;
plot(t,u(3,:));
xlabel('time');
ylabel('normal thrust');
ax = axis;
axis([t(1) t(N) -1.1 1.1]);

subplot(4,1,4);
hold on;
plot(t,norm_u);
plot(t,Lm,'r--');
xlabel('time');
ylabel('Norm of thrust');
ax = axis;
axis([t(1) t(N) -0.1 1.1]);

% trajectory
step = 1;
W = 1+ex.*cos(L)+ey.*sin(L);
Z = hx.*sin(L)-hy.*cos(L);
C = 1+hx.^2+hy.^2;
r = zeros(3,N);
r(1,:) = P.*( (1+hx.^2-hy.^2).*cos(L) + 2*hx.*hy.*sin(L) ) ./ (C.*W);
r(2,:) = P.*( (1-hx.^2+hy.^2).*sin(L) + 2*hx.*hy.*cos(L) ) ./ (C.*W);
r(3,:) = 2*P.*Z ./ (C.*W);
v(1,:) = (1./C).*sqrt(mu0./P).*( 2*hx.*hy.*(ex+cos(L))...,
         - (1+hx.^2-hy.^2).*(ey+sin(L)) );
v(2,:) = (1./C).*sqrt(mu0./P).*( (1-hx.^2+hy.^2).*(ex+cos(L))...,
         - 2*hx.*hy.*(ey+sin(L)) );
v(3,:) = (2./C).*sqrt(mu0./P).*(hx.*(ex+cos(L))+hy.*(ey+sin(L)));
Q = r;
for i = 1:N Q(:,i) = Q(:,i) / norm(Q(:,i)); end;
W = cross(r,v); 
for i = 1:N W(:,i) = W(:,i) / norm(W(:,i)); end;
S = cross(W,Q);
uc = u;
for i = 1:N
   uc(:,i) = [ Q(:,i) S(:,i) W(:,i) ] * u(:,i);
end;

% Draw Earth
rt = 6.128 ;

for i = 1:50 theta1(i) = i*2*pi/50; end;
for i = 1:50 theta2(i) = i*2*pi/50; end;

for i = 1:50
for j = 1:50
e(1,(i-1)*50+j) = rt*cos(theta1(i)).*cos(theta2(j));
e(2,(i-1)*50+j) = rt*sin(theta1(i)).*cos(theta2(j));
e(3,(i-1)*50+j) = rt*sin(theta2(j));
end;
end;

figure(fig2);
%set(fig2,'Position',[706 464 428 353]);
title(['Trajectory , T_{max} = ' num2str(Tmax) ' N, t_f = ' num2str(tf) ...
       ' h, L_f = ' num2str(L(N)) ' rad, m_f = ' num2str(M(N)) ' kg.']);
hold on;
subplot(5,1,1:3);
plot3(r(1,:),r(2,:),r(3,:));
i = 1:step:N;
hold on;
quiver3(r(1,i),r(2,i),r(3,i),uc(1,i),uc(2,i),uc(3,i),'r');
plot3(e(1,:),e(2,:),e(3,:));
hold off;
xlabel('r_1');
ylabel('r_2');
zlabel('r_3');
axis equal;
daxes3(0,0,0,'r');
drawnow;
zoom on;

subplot(5,2,[7 9]);
plot(r(1,:),r(2,:));
i = 1:step:N;
hold on;
%quiver(r(1,i),r(2,i),uc(1,i),uc(2,i));
hold off;
xlabel('r_1');
ylabel('r_2');
axis equal;
daxes(0,0,'r');
drawnow;
zoom on;

subplot(5,2,[8 10]);
plot(r(2,:),r(3,:));
i = 1:step:N;
hold on;
%quiver(r(2,i),r(3,i),uc(2,i),uc(3,i));
hold off;
xlabel('r_2');
ylabel('r_3');
daxes(0,0,'r');
axis([-47 47 -2.4 2.4]);
drawnow;
zoom on;

% State and costate
figure(fig3);
%State
subplot(7,2,1);
title(['T_{max} = ' num2str(Tmax) ' N, ctf = ' num2str(par(1)) ...
       ' m_f = ' num2str(M(N)) ' kg.']);
hold on;
xlabel ('t');
ylabel ('P');
plot(t,P);
ylabel('P');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,3);
plot(t(:),ex(:));
hold on;
plot([t(1) t(N)],[0 0],'--');
ylabel('e_x');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,5);
plot(t,ey);
hold on;
plot([t(1) t(N)],[0 0],'--');
ylabel('e_y');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,7);
plot(t,hx);
hold on;
plot([t(1) t(N)],[0 0],'--');
ylabel('h_x');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,9);
plot(t,hy);
hold on;
plot([t(1) t(N)],[0 0],'--');
ylabel('h_y');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,11);
plot(t,mod(L,2*pi));
hold on;
ylabel('L');
plot([t(1) t(N)],[pi pi],'--');
axis([t(1) t(N) 0 2*pi]);
zoom on;

subplot(7,2,13);
plot(t(:),M(:));
ylabel('m');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

% Costate
subplot(7,2,2);
hold on;
plot(t(:),p(1,:));
ylabel('p_P');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,4);
plot(t(:),p(2,:));
hold on;
plot([t(1) t(N)],[0 0],'--');
ylabel('p_{ex}');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,6);
plot(t(:),p(3,:));
hold on;
plot([t(1) t(N)],[0 0],'--');
ylabel('p_{ey}');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,8);
plot(t(:),p(4,:));
ax = axis;
hold on;
ylabel('p_{hx}');
plot([t(1) t(N)],[0 0],'--');
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,10);
plot(t(:),p(5,:));
ax = axis;
hold on;
ylabel('p_{hy}');
plot([t(1) t(N)],[0 0],'--');
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,12);
plot(t(:),p(6,:));
hold on;
ylabel('p_L');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

subplot(7,2,14);
plot(t(:),p(7,:));
ylabel('p_m');
ax = axis;
axis([t(1) t(N) ax(3) ax(4)]);
zoom on;

end;
end;
