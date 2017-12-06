%% Q1
clc; clear variables; close all;
mc=5;
T0=100;
Tw=20;
k=13;
cc=385;
odefun = @(t,Tc) k/mc/cc*(Tw-Tc); %10pt
event = @(t,Tc) ndgrid(Tc-25,1,0); %5pt
opt=odeset('event',event); %5pt
[t,Tc,tE]=ode45(odefun,[0 10*60],T0,opt); %10pt
plot(t,Tc)
fprintf('The temperature of the copper ball is 25 deg after %.1f minutes\n',tE/60) %5pt
%% Q2
syms T(t)
sol=dsolve(diff(T)==k/mc/cc*(Tw-T),T(0)==100); %10
ezplot(sol,[0 600]) %5
tSol = eval(solve(sol==25)); %5
hold on
plot(tSol,25,'o') %5
hold off
ylabel('T_c') %5
title([])
fprintf('The temperature of the copper ball is 25 deg after %.1f minutes\n',tSol/60)
%% Q3
mw=10;
cw=4179;
odefun = @(t,T) [k/mc/cc*(T(2)-T(1)); k/mw/cw*(T(1)-T(2))]; %10
[t,T]=ode45(odefun,[0 15*60],[100 20]); %10
dT = T(:,1)-T(:,2);
tSol = interp1(dT,t,1); %10
plot(t,dT,tSol,1,'o')

