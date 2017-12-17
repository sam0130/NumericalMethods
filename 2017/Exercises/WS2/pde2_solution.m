clear all
close all
clc

% parameters
L = 1500;
dx = L/400;
x = (dx/2:dx:L-dx/2)'; % ode15s requires a row vector 
l = @(x) zeros(size(x)); % no diffusion
m = @(x) ones(size(x)); % unit mass 
s = @(x) zeros(size(x)); % no source term
u = @(x) 1.5*ones(size(x)); % flow velocity
phil = 0; %concentration at x=0
phir = 0; %concentration at x=L
phi0 = 2 * (x>50 & x<250); % initial values
tmax = 600;
%% Q1
figure(1)
dt = 0.9*dx/max(u(x));
fprintf('nt %d\n',tmax/dt);
phi = phi0;
for t=dt:dt:tmax
   phi = phi+dt * convectionDiffusion(t,phi,x,m,l,s,u,phil,phir);
end
phiExact = 2 * (x-u(0)*tmax>50 & x-u(0)*tmax<250); % initial values
figure()
plot(x,phiExact,'b',x,phi,'r')
xlabel('x')
ylabel('\phi')
legend('exact','numerical')
%% Q2
figure()
phi=phi0;
for t=dt:dt:tmax
   phi = phi+dt * convectionDiffusionLW(t,phi,x,m,l,s,u,phil,phir,dt);
end
plot(x,phiExact,'b',x,phi,'r')
xlabel('x')
ylabel('\phi')
legend('exact','numerical')
%% Q4
dtConvection = 0.9*dx/max(u(x))
dtDiffusion = 0.15*dx.^2/2/max(l(x)./m(x))
%For my choice of dx, the two parameters are very similar.
%For larger dx, dtDiffusion is dominant
%% Q5
figure(3)
l = @(x) 0.5*ones(size(x));
phi=phi0;
dt = min(dtConvection,dtDiffusion);
for t=dt:dt:tmax
   phi = phi+dt * convectionDiffusion(t,phi,x,m,l,s,u,phil,phir);
end
figure()
plot(x,phi0,x,phi,'r')
xlabel('x')
ylabel('\phi')
legend('t=0',num2str(tmax,'t=%ds'))
legend('initial concentration','final concentration')
%% Q6
W = 50;
D = 4;
T0 = trapz(x,phi0)*W*D
T = trapz(x,phi)*W*D
error=(T0-T)/T
disp('The total concentration is almost perfectly conserved')
%% Q7
W = 100;
dy =W/10;
y = (dy/2:dy:W-dy/2)';
% create a mesh
[y,x] = meshgrid(y,x);
% define all functions in 2D
l = @(x,y) 0.5*ones(size(x));
m = @(x,y) ones(size(x));
s = @(x,y) zeros(size(x));
u = @(x,y) 0.8+1.2*((W-y).*y)/W^2;
v = @(x,y) zeros(size(x));
phi0 = 2 * (x>50 & x<250); % initial values
tmax = 600;

dt = min([0.9*dx/sqrt(max(max(u(x,y).^2+v(x,y).^2))),...
          0.5*min(dx,dy)/2/max(max(l(x,y)./m(x,y)))]);
fprintf('nt %d\n',tmax/dt);
phi2D=phi0;
for t=dt:dt:tmax
   phi2D = phi2D+dt * convectionDiffusion2DPipe(t,phi2D,x,y,m,l,s,u,v,phil);
end
figure()
surf(x,y,phi2D)
xlabel('x')
ylabel('y')
axis tight
