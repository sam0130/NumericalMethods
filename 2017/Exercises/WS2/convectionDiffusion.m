function dphi = convectionDiffusion(t,phi,x,m,l,s,u,phil,phir)
% Compute time derivative for transport equation
%     dphi/dt + u * dphi/dx + v * dphi/dy 
%             = \nabla.(\lambda\nabla\phi) + s
% For dirichlet boundary conditions
%
% input variables:
%   t - current time (scalar)
%   phi - current value (nx1 vector)
%   x - position (nx1 vector)
%   m - mass (function handle @(x,y))
%   l - diffusion coefficient (function handle @(x,y))
%   s - source (function handle @(x,y))
%   u - x-velocity (function handle @(x,y))
%   phil - b.c. (scalar)
%   phir - b.c. (scalar)
%
% output variables:
%   dphi - derivative of temperature (nx1 vector)
dx = x(2)-x(1);
% Define time derivative dT/dt
dphi = zeros(size(phi));
% Compute diffusive term by central differencing
for i=2:length(phi)-1
    dphi(i) = ( l(x(i)+dx/2)*(phi(i+1)-phi(i)) ...
            - l(x(i)-dx/2)*(phi(i)-phi(i-1)) )/dx^2 ...
            - m(x(i))*u(x(i))*(phi(i)-phi(i-1))/dx;
end
% Boundary conditions
dphi(1) = (   l(x(1)+dx/2)*(phi(2)-phi(1)) ...
        - 2*l(x(1)-dx/2)*(phi(1)-phil  ) )/dx^2 ...
        - m(x(1))*u(x(1))*(phi(1)-phil)/dx;
dphi(end) = ( 2*l(x(end)+dx/2)*(phir    -phi(end)  ) ...
            - l(x(end)-dx/2)*(phi(end)-phi(end-1)) )/dx^2 ...
            - m(x(end))*u(x(end))*(phi(end)-phi(end-1))/dx;
% Add source term and divide by mass
dphi = (dphi + s(x))./m(x);
end