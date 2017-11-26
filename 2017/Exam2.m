%%%% EXAM2 %%%%%%%%%%
clc; close all; clear variables;
%% Question1
d = 0.05;
ho = 0.5;
g = 9.8;
syms h(t)

h = dsolve( diff(h) == (-d^2/(4*h^2))*sqrt(2*g*h), h(0) == ho);
tf = eval(solve(h==0,t));
figure
ezplot(h,[0 tf])
%% Question2
D = sqrt(1/3);
ho = 0.5;
d = 0.05;

syms h(t)
h = dsolve( diff(h) == (-d^2/(D^2))*sqrt(2*g*h), h(0) == ho);
tf = eval(solve(h(2)==0,t));
figure
ezplot(h(2),[0 tf(2)])
%% Question 3
d = 2;
To = 100;
Ts = 20;
k = 0.004;
tf = 100;
dt = 0.01;
odefun = @(t,T) -k*(T-Ts);
[t,y] = odeEuler(odefun, [0 tf], To, dt);
Tf = y(end);
fprintf('The temperature after tf=%.2f s is Tf=%.2f \n', tf, Tf);

%% Question 4
m = 2;
k = 10;
gamma = 1;
y0 = [0; 0.5];
T = [0 20];
dt = 0.1;

event = @(t,y) ndgrid(y(2)-0, 0, 0 );
opt = odeset('Event',event);

odeSpring = @(t,y) [y(2); (-k/m)*y(1) - (gamma/m)*y(2)];
[t,y,te,ye] = ode45(odeSpring, T,y0, opt);
plot(t, y(:,1), te, ye(:,1), 'o')

%% Question 5
mu = 100;
x0 = [1; 0];
T = [0 200];

odeVanDerPol = @(t,x) [x(2); -x(1) - mu*(x(1)^2 - 1)*x(2)];
[t_ode45,y_ode45] = ode45(odeVanDerPol, T,x0);
[t_ode15s,y_ode15s] = ode15s(odeVanDerPol, T,x0);

figure
hold on
plot(t_ode45, y_ode45(:,1))
plot(t_ode15s, y_ode15s(:,1))
legend('ode45','ode15s')

