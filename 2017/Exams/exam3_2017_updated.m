%% Wall
clear variables
close all
clc
l2 = 1.2; l1 = 0.5;
l = @(x) l1+(l2-l1)*(x>=.05);
rc1 = 500*2000; rc2 = 1500*1000;
rc = @(x) rc1+(rc2-rc1)*(x>=.05);
L = .30;
tMax = 86400*5;
To = @(t) 15+10*sin(2*pi*t/86400);
Ti = 25;

pdefun = @(x,t,T,dTdx) ndgrid(rc(x),l(x)*dTdx,0);
icfun = @(x) Ti;
bcfun = @(xl,Tl,xr,Tr,t) ndgrid(Tl-To(t),0,Tr-Ti,0);
xmesh = linspace(0,L);
tmesh = linspace(0,tMax, 500);
T=pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);

figure()
plot(xmesh,T(1:50:end,:))
xlabel('x')
ylabel('T')
axis tight

for i=1:length(tmesh)
    phi0=gradient(T(i,:),xmesh);
    phi(i)=phi0(end);
end
figure()
plot(tmesh/86400,phi)
xlabel('t')
ylabel('\phi(x=0.3)')
axis tight

%% Couette Cell
clc
clear variables
close all
mu = 1e-3; %Pa.s = kg/m/s (nu=mu/rho = 3e-6 m^2/s)
rho = 1e3; %kg/m^3 (1 g/m^3)
R1 = 9e-3;
R2 = 10e-3;
tMax = 10;
wo=0.1;

icfun = @(x) 0;
bcfun = @(xl,wl,xr,wr,t) ndgrid(wl,0,wr-wo,0);
xmesh = linspace(R1,R2);
tmesh = linspace(0,tMax,50);
pdefun = @(x,t,w,dwdx) ndgrid(rho,mu*dwdx,0);
w=pdepe(1,pdefun,icfun,bcfun,xmesh,tmesh);
figure()
plot(xmesh,w(end,:))
axis tight
xlabel('r')
ylabel('v')
tau=mu*xmesh.*gradient(w(end,:),xmesh);
%fprintf('The shear stress on the wall is %f ')

figure()
plot(xmesh,tau)

%% Tracer in a pipe
clc
clear variables
close all
D = 0.1;
Diff = 0.0238;
w = 1;
L = 400;

dx = L/1000;

x = 0:dx:L;

l = @(x) Diff*ones(size(x)); 
m = @(x) ones(size(x)); % unit mass 
s = @(x) zeros(size(x)); % no source term
u = @(x) 1*ones(size(x)); % flow velocity


Cl = 0; %concentration at x=0
Cr = 0; %concentration at x=L
Co = @(x) (0.001/sqrt(2*pi*w))*exp(-(x.^2)/2*w^2); % initial

tmax = 5*60;

dtConvection = 0.9*dx/max(u(x));
dtDiffusion = 0.02*dx.^2/2/max(l(x)./m(x));

C = Co(x);
dt = min(dtConvection,dtDiffusion);
for t=dt:dt:tmax
   C = C + dt * convectionDiffusion(t,C,x,m,l,s,u,Cl,Cr);
end
plot(x,C)
%% analytical
clc;
close all;
A = pi*D^2/4;
x = 0:dx:L;
t = 0:dt:tmax;

[T,X] = meshgrid(x,t);
C_anal = @(x,t) w/A/1000./sqrt(4*pi*Diff*t).*exp(- (x-u(x).*t).^2./(4*Diff*t));
figure()
surf(x,t,C_anal(X,T))
xlabel('x')
ylabel('t')
figure()
hold on
plot(x,C)
plot(x,C_anal(x,tmax))





