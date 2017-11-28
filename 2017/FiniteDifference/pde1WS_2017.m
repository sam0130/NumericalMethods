clear variables;
clc;
close all;
%% Q1 State the IBVP
% Find T(t,x) on 0<t<tf and 0<x<L such that
% d/dt(r*C*T) = d/dx(l*d/dTx)
% with T(0,x)=To and T(t,0)=To, T(t,L)=Ti
%% Q2 Finite difference formulation
% the implementation of dT is at the end of the file

% density, thermal capacity, heat transfer coefficient 
% (function handles @(x))
x1 = 0.1;
rb = 1500; cb = 1000; lb = 1.2;
ri = 500; ci = 2000; li = 0.1;
r = @(x) ri + (rb-ri)*heaviside(x-x1);          % density
c = @(x) ci + (cb-ci)*heaviside(x-x1);          % Cp
l = @(x) li + (lb-li)*heaviside(x-x1);          % K conductivity
% mass, source (function handles @(x))
m = @(x) r(x).*c(x);                            % rho*Cp
s = @(x) zeros(size(x));                        
% time interval
tmax = 3600*24*4;
% spatial discretisation
x = linspace(0,.4,40)';                        
x = x(1:end-1) + diff(x)/2;                    % center points  
% initial conditions (nx1 vector)
T0 = -10*ones(size(x));
% boundary conditions (scalar)
Tl = -10;
Tr = 25;
% Time step based on diffusion condition
dx = x(2)-x(1);
a = max(l(x)./m(x));
dt = 0.5*dx*dx/2/a;
%% Q3 Forward Euler solution
T = T0;
for t = 0:dt:tmax
    T=T+dt*odefun(t,T,x,m,l,s,Tl,Tr);
end
plot(x,T);
xlabel('x'); ylabel('T');
%% Q4 local heat flux
phi = l(x).*gradient(T,x);
plot(x,phi,'.')
xlabel('x'); ylabel('\phi');
%% Q5 Heat loss
fprintf('Heat loss is %.0f Watt\n',median(phi)*6*2.5)
%% Q7 ode15s solution
odefun2 = @(t,T) odefun(t,T,x,m,l,s,Tl,Tr);
[t15,T15] = ode15s(odefun2,linspace(0,tmax,5),T0);
plot(x,T,x,T15(end,:),'.')
legend('Euler','ode15s');
xlabel('x'); ylabel('T');
%% Q8 ode15s vs ode45
tic
[t,T] = ode45(odefun2,linspace(0,tmax,5),T0);
toc
tic
[t,T] = ode15s(odefun2,linspace(0,tmax,5),T0);
toc
%% Q9
% mass, source (function handles @(x))
m = @(x) 8960*385*ones(size(x));
l = @(x) 400*ones(size(x));
s = @(x) 80000 * heaviside(x-0.4).*heaviside(0.6-x);
% time interval
tmax = 3600*1;
% boundary conditions (scalar)
Tl = 20;
Tr = 20;
% spatial discretisation
x = linspace(0,1,50)';
x = x(1:end-1) + diff(x)/2;
% initial conditions (nx1 vector)
T0 = 20*ones(size(x));
odefun2 = @(t,T) odefun(t,T,x,m,l,s,Tl,Tr);
[t,T] = ode15s(odefun2,linspace(0,tmax,7),T0);
plot(x,T)
legend(num2str(t/60,'t=%5.0fm'));
xlabel('x'); ylabel('T');
%%
odefun2 = @(t,T) odefunNeumann(t,T,x,m,l,s,Tr);
[t,T] = ode15s(odefun2,linspace(0,3600*6,7),T0);
plot(x,T)
legend(num2str(t/3600,'t=%5.0fh'));
xlabel('x'); ylabel('T');

