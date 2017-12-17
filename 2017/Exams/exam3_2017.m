%% Wall
clear variables
clc
l1 = 1.2; l2 = 0.5;
l = @(x) l1+(l2-l1)*(x>=.05);
rc1 = 1500*1000; rc2 = 500*2000;
rc = @(x) rc1-(rc2-rc1)*(x>=.05);
L = .25;
tMax = 86400*5;
To = @(t) 15-5*cos(2*pi*t/86400);
Ti = 25;

pdefun = @(x,t,T,dTdx) ndgrid(rc(x),l(x)*dTdx,0);
icfun = @(x) 25;
bcfun = @(xl,Tl,xr,Tr,t) ndgrid(Tl-To(t),0,Tr-Ti,0);
xmesh = linspace(0,L);
tmesh = linspace(0,tMax);
T=pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);

figure(1)
plot(xmesh,T([90 100],:))
xlabel('x')
ylabel('T')
legend('12:00','24:00')
axis tight
figure(2)
for i=1:length(tmesh)
    phi0=gradient(T(i,:),xmesh);
    phi(i)=phi0(end);
end
plot(tmesh,phi)
xlabel('t')
ylabel('\phi(x=0.25)')
axis tight


%% Couette cell
clc
clear variables
close all
mu = 1e-3; %Pa.s = kg/m/s (nu=mu/rho = 3e-6 m^2/s)
rho = 1e3; %kg/m^3 (1 g/m^3)
R1 = 10e-3;
R2 = 12e-3;
tMax = 10;
wo=0.1;

icfun = @(x) 0;
bcfun = @(xl,wl,xr,wr,t) ndgrid(wl,0,wr-wo,0);
xmesh = linspace(R1,R2);
tmesh = linspace(0,tMax, 50);
pdefun = @(x,t,v,dwdx) ndgrid(rho,mu*dwdx,0);
w=pdepe(1,pdefun,icfun,bcfun,xmesh,tmesh);
figure(3)
plot(xmesh,w(end,:).*xmesh)
axis tight
xlabel('r')
ylabel('v')
tau=mu*xmesh.*gradient(w(end,:),xmesh);
fprintf('The shear stress on the wall is %f ')
plot(tau)

%% Sphere
clear all
clc
l = 400;
rc = 8960*385;
tMax = 12;
Tw = 20;
Tc = 100;
R = 0.06;

pdefun = @(x,t,T,dTdx) ndgrid(rc,l*dTdx,0);
icfun = @(x) Tc;
bcfun = @(xl,Tl,xr,Tr,t) ndgrid(0,1,Tr-Tw,0);
xmesh = linspace(0,R);
tmesh = linspace(0,tMax);
T=pdepe(2,pdefun,icfun,bcfun,xmesh,tmesh);

figure(4)
plot(tmesh,T(:,1))
xlabel('x')
ylabel('T')
axis tight

t30 = interp(T(:,1),tmesh,30)
