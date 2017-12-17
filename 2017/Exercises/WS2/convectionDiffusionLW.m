function dphi = convectionDiffusion(t,phi,x,m,lorig,s,v,phil,phir,dt)
% Compute time derivative dT/dt for static heat equation
%     dT/dt = d/dx(\lambda dT/dx) + s
% input variables:
%   t - current time (scalar)
%   T - current temperature (nx1 vector)
%   x - position (nx1 vector)
%   m - mass (function handle @(x))
%   l - diffusion coefficient (function handle @(x))
%   s - source (function handle @(x))
%   Tl - left b.c. (scalar)
%   Tr - right b.c. (scalar)
%
% output variables:
%   dT - derivative of temperature (nx1 vector)
dx = x(2)-x(1);
% Define time derivative dT/dt
dphi = zeros(size(phi));
% Compute diffusive term by central differencing
l=@(x) lorig(x)+v(x).^2*dt/2;
for i=2:length(phi)-1
    dphi(i) = ( l(x(i)+dx/2)*(phi(i+1)-phi(i)) ...
            - l(x(i)-dx/2)*(phi(i)-phi(i-1)) )/dx^2 ...
            - v(x(i))*(phi(i+1)-phi(i-1))/dx/2;
end
% Boundary conditions
dphi(1) = (   l(x(1)+dx/2)*(phi(2)-phi(1)) ...
        - 2*l(x(1)-dx/2)*(phi(1)-phil  ) )/dx^2 ...
        - v(x(1))*(phi(2)-phil)/dx/1.5;
dphi(end) = ( 2*l(x(end)+dx/2)*(phir    -phi(end)  ) ...
             - l(x(end)-dx/2)*(phi(end)-phi(end-1)) )/dx^2 ...
             - v(x(end))*(phir-phi(end-1))/dx/1.5;
% Add source term and divide by mass
dphi = (dphi + s(x))./m(x);
end