%% Flow in a blood vessel
close all
clear variables
clc

%% Q1 constant blood flow
%pg 149, glycerol-water mixture, close to blood
mu = 12e-3; %Pa.s = kg/m/s (nu=mu/rho = 3e-6 m^2/s)
rho = 1e3; %kg/m^3 (1 g/m^3)
A = 5000; %15*133; % N/m (1 dyne/cm)
R = 2.7e-3;
tMax = 100;

% format: pdefun = @(x,t,v,dvdx) ndgrid(c,f,s);
% equation: d/dt(rho v) = d/dy(mu*dv/dy) - dp/dx;
pdefun = @(x,t,v,dvdx) ndgrid(rho,mu*dvdx,A);

% format: icfun = @(x) v0;
icfun = @(x) 0;

% format: bcfun = @(xl,vl,xr,vr,t) ndgrid(pl,ql,pr,qr)
% no-slip: v(xl)=v(xb)=0
bcfun = @(xl,vl,xr,vr,t) ndgrid(vl-0,0,vr-0,0);

xmesh = linspace(0,R);

tmesh = linspace(0,tMax);

opt = odeset('AbsTol',1e-10);
v=pdepe(1,pdefun,icfun,bcfun,xmesh,tmesh,opt);
surf(xmesh,tmesh,v)
view(90,0)
xlabel('x')
ylabel('t')
zlabel('v')
plot([-xmesh(end:-1:1) xmesh]*100,[v(end,end:-1:1) v(end,:)]*100)
xlabel('x [cm]')
ylabel('v_z [cm/s]')
axis tight
title('steady state velocity profile')
set(gcf,'Position',[900 0 300 200])
saveas(gcf,'blood1.png')

%% Q2 measure the Reynolds number
vmax = max(v(end,:));
Re = rho*vmax*R/mu

%% Q3 measure blood flow
Qmax = trapz(xmesh,v(end,:)*2*pi.*xmesh)

%% Q4 oscillating blood flow
w = 6*pi; % 3 Hz
tmesh = linspace(0,1.1,111);
pdefun = @(x,t,v,dvdx) ndgrid(rho,mu*dvdx,A*cos(w*t));
v=pdepe(1,pdefun,icfun,bcfun,xmesh,tmesh,opt);
% axial velocity profile
plot(tmesh,v(:,1)*100)
xlabel('t [s]')
ylabel('v_z(r=0) [cm/s]')
title('axial velocity profile')
xlim([0,1.1])
set(gcf,'Position',[900 0 300 200])
saveas(gcf,'blood2.png')
%% Q5 radial pressure profile
id=91:2:107
figure()
plot([-xmesh(end:-1:1) xmesh]*100,[v(id,end:-1:1) v(id,:)]*100)
xlabel('r [cm]')
ylabel('v_z [cm/s]')
legend(num2str(tmesh(id)','%.2fs'))
axis tight
title('radial velocity profile')
set(gcf,'Position',[900 0 300 200])
saveas(gcf,'blood3.png')

%% Q6 oscillating blood flow, no back flow
pdefun = @(x,t,v,dvdx) ndgrid(rho,mu*dvdx,max(0,A*cos(w*t)));
v=pdepe(1,pdefun,icfun,bcfun,xmesh,tmesh,opt);
%% axial velocity profile
plot(tmesh,v(:,1)*100)
xlabel('t [s]')
ylabel('v_z(r=0) [cm/s]')
title('axial velocity profile')
set(gcf,'Position',[900 0 300 200])
saveas(gcf,'blood4.png')
%% measure blood flow
for i = 1:length(tmesh)
    Q(i) = trapz(xmesh,v(i,:)*2*pi.*xmesh);
end
plot(tmesh,Q/Qmax)
xlabel('t')
ylabel('Q/Q_{max}')
title('mass flow rate')
set(gcf,'Position',[900 0 300 200])
saveas(gcf,'blood5.png')


