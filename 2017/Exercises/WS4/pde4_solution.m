% heat transfer coeff, bulk density x therm.cap, source
x1 = 0.2;   % Bricks
x2 = 0.31;  % brick + insulation
l1 = 1.2;   % brick heat transfer coeff
l2 = 0.1;   % insulation heat transfer coeff
l3 = 0.7;   % concrete heat transfer
l = @(x) l1+(l2-l1)*heaviside(x-x1)+(l3-l2)*heaviside(x-x2);
c1 = 1500*1000; c2 = 500*2000; c3 = 2000*1000;  % rho*Cp                  
c = @(x) c1+(c2-c1)*heaviside(x-x1)+(c3-c2)*heaviside(x-x2);
s = 0;

% define pde function c*dTdt = d/dx(f)+s (for m=0):
pdefun = @(x,t,T,dTdx) ndgrid(c(x),l(x)*dTdx,s);

% initial conditions for temperature
icfun = @(x) 20;

% boundary conditions (p+q*phi=0).
bcfun = @(xl,Tl,xr,Tr,t) ndgrid(Tl-1200, 0, Tr-20, 0);

% spatial/temporal discretisation
xmesh = linspace(0,0.61);
tmesh = linspace(0,5*24*3600,6);

%solve initial-boundary value problem
T = pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);
figure()
plot(xmesh,T,'.-')
legend(num2str(tmesh'/24/3600,'t=%.0f d'))

%% Q2

%ode function
odefun=@(x,y)[y(2)/l(x); -s]; 

% boundary conditions, given as p+q*f=0.
bcfun=@(yl,yr)[yl(1)-1200;yr(1)-20]; 

% spatial/temporal discretisation
solinit=bvpinit(xmesh,[0 0]); 

% solve the initial-boundary value problem
sol=bvp4c(odefun,bcfun,solinit);
figure()
plot(sol.x,sol.y(1,:),'-',xmesh,T(end,:),'-')                              % Compare Temperature from both solvers

%% Q3
f = l(xmesh).*gradient(T(end,:),xmesh);
figure()
plot(sol.x,sol.y(2,:),'-',xmesh,f,'-')

%% Q4
Tmax = interp1(sol.x,sol.y(1,:),0.31);
fprintf('The maximum temperature in the concrete layer is %.0f deg\n',Tmax)

%% Q5
% heat transfer coeff, bulk density x therm.cap, source
x1 = 0.7; x2 = 0.81;
l1 = 1.2; l2 = 0.1; l3 = 0.7;
l = @(x) l1+(l2-l1)*heaviside(x-x1)+(l3-l2)*heaviside(x-x2);
c1 = 1500*1000; c2 = 500*2000; c3 = 2000*1000;
c = @(x) c1+(c2-c1)*heaviside(x-x1)+(c3-c2)*heaviside(x-x2);
s = 0;

% define pde function c*dTdt = d/dx(f)+s (for m=0):
pdefun = @(x,t,T,dTdx) ndgrid(c(x),l(x)*dTdx,s);

% initial conditions for temperature
icfun = @(x) 20;

% boundary conditions (p+q*phi=0).
bcfun = @(xl,Tl,xr,Tr,t) ndgrid(Tl-1200, 0, Tr-20, 0);

% spatial/temporal discretisation
xmesh = linspace(0.5,1.11);
tmesh = linspace(0,5*24*3600,6);

%solve initial-boundary value problem
T = pdepe(0,pdefun,icfun,bcfun,xmesh,tmesh);
TC = pdepe(1,pdefun,icfun,bcfun,xmesh,tmesh);
figure()
plot(xmesh,TC(end,:),'-',xmesh,T(end,:),'-')

%% Q6
Tmax = interp1(xmesh,T(end,:),0.81)
TmaxC = interp1(xmesh,TC(end,:),0.81)

