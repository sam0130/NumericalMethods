%% Clean up and define parameters
syms clear
clear all
close all
format compact
% Define parameters
C=0.44;    
rhoF=1000; 
rhoS=7900;  
r=4e-2;
gr=9.8;
g=(rhoS-rhoF)/rhoF*gr;
gam=3*rhoF*C/8/rhoS/r;

%% Question 1
% function: v(t)
% ic: v(0)=0 (5pt)

%% Question 2
% Solve IVP analytically (10pt)
syms v(t)
vSol = dsolve( diff(v) == -gam*v^2+g, v(0) == 0);
%eval(limit(vSol,t,inf))
% Solve v(t)=11 for t (5pt)
tf=eval(solve(vSol==11,t));
disp(['Question 2: tf=' num2str(tf,8) ' s'])

%% Question 3
% define odefun (5pt)
dv = @(t,v) -gam*v^2+g; 
% set relative accuracy (5pt)
opt=odeset('RelTol',1e-4);
% use ode45 to solve IVP (5pt)
[t,v]=ode45(dv,[0 tf],0,opt);
disp(['Question 3: vf=' num2str(v(end),8) ' m/s'])

%% Question 4
%use a stopevent (10pt)
opt=odeset('Event',@stopevent);
[t,v]=ode45(dv,[0 100],0,opt);
disp(['Question 4: tf=' num2str(t(end),8) ' s']);

%% Question 5
% plot both graphs (10pt)
ezplot(vSol,[0 tf])
hold all
plot(t,v,'.')
hold off
xlabel('time [s]')
ylabel('velocity [m/s]')
title('Question 5')

%% Question 6
% use forward euler (10pt)
[t,v1]=odeEuler(dv,[0 tf],0,tf/100);
disp(['Question 6: vf=' num2str(v(end),8) ' m/s'])

%% Question 6 alternative answer: self-written euler
dt=tf/100;
t=(0:dt:tf)';
v = zeros(size(t));
v(1) = 0;
for i=2:length(t)
    v(i) = v(i-1)+dt*(-gam*v(i-1)^2+g);
end
disp(['Question 6: vf=' num2str(v(end),8) ' m/s']);

%% Question 7
% function: y(t), v(t) (5pt)
% ic: y(0)=0, v(0)=0 (5pt)

%% Question 8
% define odefun (5pt)
dy = @(t,y) [y(2); -gam*y(2)^2+g]; 
% set relative accuracy (5pt)
opt=odeset('RelTol',1e-4);
% use ode45 to solve IVP (5pt)
[t,y]=ode45(dy,[0 tf],[0 0],opt);
z=y(:,1);
disp(['Question 8: zf=' num2str(z(end)) ' m']);


