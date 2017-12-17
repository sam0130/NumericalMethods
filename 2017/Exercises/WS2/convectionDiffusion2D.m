function dphi = convectionDiffusion2D(t,phi,x,y,m,l,s,u,v,phib)
% Compute time derivative for transport equation
%     dphi/dt + u * dphi/dx + v * dphi/dy 
%             = \nabla.(\lambda\nabla\phi) + s
% for Dirichlet boundary conditions.
%
% input variables:
%   t - current time (scalar)
%   phi - current value (nxm matrix)
%   x - position (nxm matrix)
%   y - position (nxm matrix)
%   m - mass (function handle @(x,y))
%   l - diffusion coefficient (function handle @(x,y))
%   s - source (function handle @(x,y))
%   u - x-velocity (function handle @(x,y))
%   v - y-velocity (function handle @(x,y))
%   phib - b.c. (function handle @(x,y))
%
% output variables:
%   dphi - derivative of temperature (nxm matrix)
dx = x(2,1)-x(1,1);
dy = y(1,2)-y(1,1);
% Define time derivative dT/dt
dphi = zeros(size(phi));
% Compute diffusive term by central differencing
for i=2:size(phi,1)-1
    for j=2:size(phi,2)-1
        dphi(i,j) = ( l(x(i,j)+dx/2,y(i,j))*(phi(i+1,j)-phi(i,j)) ...
                    - l(x(i,j)-dx/2,y(i,j))*(phi(i,j)-phi(i-1,j)) )/dx^2 ...
                  + ( l(x(i,j),y(i,j)+dy/2)*(phi(i,j+1)-phi(i,j)) ...
                    - l(x(i,j),y(i,j)-dy/2)*(phi(i,j)-phi(i,j-1)) )/dy^2 ...
                   - m(x(i,j),y(i,j))*u(x(i,j),y(i,j))*(phi(i,j)-phi(i-1,j))  /dx ...
                   - m(x(i,j),y(i,j))*v(x(i,j),y(i,j))*(phi(i,j)-phi(i,j-1))  /dy;
    end
end
for j=2:size(phi,2)-1
    i=1;
    dphi(i,j) = ( l(x(i,j)+dx/2,y(i,j))*(phi(i+1,j)-phi(i,j)) ...
                - l(x(i,j)-dx/2,y(i,j))*(phi(i,j)-phib) )/dx^2 ...
              + ( l(x(i,j),y(i,j)+dy/2)*(phi(i,j+1)-phi(i,j)) ...
                - l(x(i,j),y(i,j)-dy/2)*(phi(i,j)-phi(i,j-1)) )/dy^2 ...
               - m(x(i,j),y(i,j))*u(x(i,j),y(i,j))*(phi(i,j)-phib)        /dx ...
               - m(x(i,j),y(i,j))*v(x(i,j),y(i,j))*(phi(i,j)-phi(i,j-1))  /dy;
    i=size(phi,1);
    dphi(i,j) = ( l(x(i,j)+dx/2,y(i,j))*(phib-phi(i,j)) ...
                - l(x(i,j)-dx/2,y(i,j))*(phi(i,j)-phi(i-1,j)) )/dx^2 ...
              + ( l(x(i,j),y(i,j)+dy/2)*(phi(i,j+1)-phi(i,j)) ...
                - l(x(i,j),y(i,j)-dy/2)*(phi(i,j)-phi(i,j-1)) )/dy^2 ...
                - m(x(i,j),y(i,j))*u(x(i,j),y(i,j))*(phi(i,j)-phi(i-1,j))  /dx ...
                - m(x(i,j),y(i,j))*v(x(i,j),y(i,j))*(phi(i,j)-phi(i,j-1))  /dy;
end
for i=2:size(phi,1)-1
    j=1;
    dphi(i,j) = ( l(x(i,j)+dx/2,y(i,j))*(phi(i+1,j)-phi(i,j)) ...
                - l(x(i,j)-dx/2,y(i,j))*(phi(i,j)-phi(i-1,j)) )/dx^2 ...
              + ( l(x(i,j),y(i,j)+dy/2)*(phi(i,j+1)-phi(i,j)) ...
                - l(x(i,j),y(i,j)-dy/2)*(phi(i,j)-phib) )/dy^2 ...
                - m(x(i,j),y(i,j))*u(x(i,j),y(i,j))*(phi(i,j)-phi(i-1,j))  /dx ...
                - m(x(i,j),y(i,j))*v(x(i,j),y(i,j))*(phi(i,j)-phib)        /dy;
    j=size(phi,2);
    dphi(i,j) = ( l(x(i,j)+dx/2,y(i,j))*(phi(i+1,j)-phi(i,j)) ...
                - l(x(i,j)-dx/2,y(i,j))*(phi(i,j)-phi(i-1,j)) )/dx^2 ...
              + ( l(x(i,j),y(i,j)+dy/2)*(phib-phi(i,j)) ...
                - l(x(i,j),y(i,j)-dy/2)*(phi(i,j)-phi(i,j-1)) )/dy^2 ...
                - m(x(i,j),y(i,j))*u(x(i,j),y(i,j))*(phi(i,j)-phi(i-1,j))  /dx ...
                - m(x(i,j),y(i,j))*v(x(i,j),y(i,j))*(phi(i,j)-phi(i,j-1))  /dy;
end
% Add source term and divide by mass
dphi = (dphi + s(x))./m(x);
end