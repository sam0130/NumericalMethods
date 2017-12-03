%% Question 1: recreate fig. 7.8. using the properties of concrete
clc; clear variables; close all;

L = 0.1; % wall thickness
xmesh = linspace(-L,L); % spatial discretisation
tmax = 3600;
tmesh = linspace(0,tmax,7); % temporal discretisation

lambda = 1.35; % heat transfer coefficient
rho = 1000; % density
C = 2000; % thermal capacity
q = 0; % heat source
       
% define pde function c*dTdt = d/dx(f)+s (for m=0):
% c bulk density x therm.cap.
% f heat flux
% s heat source 
pdefun = @(x,t,T,dTdx) ndgrid(rho*C, lambda*dTdx, q);

T0 = 20;
% initial conditions for temperature
icfun = @(x) T0;

TB = 30; %Temperature of 20 degree inside the house
% boundary conditions (p+q*phi=0).
% For Dirichlet b.c. T=TB, set p = T-Tb, q = 0
bcfun = @(xl,Tl,xr,Tr,t) ndgrid(Tl-TB, 0, Tr-TB, 0);

%solve the initial-boundary value problem
T = pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);
figure()
plot(xmesh,T(end,:),'.-')

%% Question 2
figure()
plot(xmesh,T,'.-')
xlabel('x')
ylabel('T')
legend(num2str(tmesh'/60,'t=%.f m')) %nice way to create a legend
axis tight

%% Question 3: repeat using Neumann bc at x=0

xmesh = linspace(0,L); % new spatial discretisation

% new boundary conditions (p+q*phi=0).
% For Neumann b.c. phi=0, set p=0, q=1
bcfun = @(xl,Tl,xr,Tr,t) ndgrid(0, 1, Tr-TB, 0); 

%solve the initial-boundary value problem
T = pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);

figure()
plot(xmesh,T,'.-')
xlabel('x')
ylabel('T')
legend(num2str(tmesh'/3600,'t=%.0f h'))
axis tight

%% Question 3

% solve until F0=2
F0max=2;
a=lambda/rho/C;
tmax=F0max/a*L^2;

% write more time steps
tmesh = linspace(0,tmax,50);

% solve
T = pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);

% compute mean and center values
TM = mean(T,2);
TG = T(:,1);

% to test, plot the mean and center values
figure()
plot(tmesh,TM,'.-',tmesh,TG,'.-');
xlabel('t')
legend('T_M','T_G')

%% plot the non-dimensionalised values
F0=a*tmesh/L^2;
YM=(TB-TM)./(TB-T0);
YG=(TB-TG)./(TB-T0);
figure()
semilogy(F0,YM,'.-')
hold on
semilogy(F0,YG,'.-');
xlabel('F_0')
ylabel('Y')
legend('Y_M','Y_G')

%% Question 4
    
%Add the first term in the Fourier series to the plot()

hold on
YF=2/(0.5*pi)*exp(-(0.5*pi)^2*F0)*mean(cos(0.5*pi*xmesh/L))...
   -2/(1.5*pi)*exp(-(1.5*pi)^2*F0)*mean(cos(1.5*pi*xmesh/L));
semilogy(F0,YF,'-x')
legend('Y_M','Y_G','Y_G^F')
hold off

%% Question 5

% go back to original time interval
tmax = 3600;
tmesh = linspace(0,tmax,5);

% solve for different geometries
TSlab = pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);
TCylinder = pdepe(1,pdefun,icfun,bcfun,xmesh,tmesh);
TSphere = pdepe(2,pdefun,icfun,bcfun,xmesh,tmesh);
plot(xmesh,TSlab(end,:),xmesh,TCylinder(end,:),xmesh,TSphere(end,:));
xlabel('x')
ylabel('T')
legend('Slab','Cylinder','Sphere')

%% Question 6

% solve until F0=2
F0max=2;
a=lambda/rho/C;
tmax=F0max/a*L^2;

% write more time steps
tmesh = linspace(0,tmax,50);

% solve
TSlab = pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);
opt = odeset('RelTol',1e-8);
TCylinder = pdepe(1,pdefun,icfun,bcfun,xmesh,tmesh,opt);

% compute mean and center values
TMSlab = mean(TSlab,2);
TGSlab = TSlab(:,1);
for i=1:size(TCylinder,1)
    TMCylinder(i) = trapz(xmesh,xmesh.*TCylinder(i,:))./trapz(xmesh,xmesh);
end
TGCylinder = TCylinder(:,1);

% to test, plot the mean and center values
figure()
plot(tmesh,TMSlab,'.-',tmesh,TGSlab,'.-',...
    tmesh,TMCylinder,'x-',tmesh,TGCylinder,'x-');
xlabel('t')
legend('T_M^{slab}','T_G^{slab}','T_M^{cyl}','T_G^{cyl}')

%% plot the non-dimensionalised values
F0=a*tmesh/L^2;
YM=(TB-TM)./(TB-T0);
YG=(TB-TG)./(TB-T0);
YMCylinder=(TB-TMCylinder)./(TB-T0);
YGCylinder=(TB-TGCylinder)./(TB-T0);
semilogy(F0,YM,'.-',F0,YG,'.-',F0,YMCylinder,'.-',F0,YGCylinder,'.-');
xlabel('F_0')
ylabel('Y')
legend('Y_M^{slab}','Y_G^{slab}','Y_M^{cyl}','Y_G^{cyl}')

