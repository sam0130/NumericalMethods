function dT = odefunNeumann(t,T,x,m,l,s,Tr)
% computes temperature derivative based on heat equation
%
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
dT = zeros(size(T));

% Compute the diffusive term using central differencing
for i=2:length(T)-1
    dT(i) = ( l(x(i)+dx/2)*(T(i+1)-T(i)) ...
            - l(x(i)-dx/2)*(T(i)-T(i-1)) )/dx^2;
end

% Boundary conditions
dT(1) = (   l(x(1)+dx/2)*(T(2)-T(1)))/dx^2;
dT(end) = ( 2*l(x(end)+dx/2)*(Tr    -T(end)  ) ...
            - l(x(end)-dx/2)*(T(end)-T(end-1)) )/dx^2;

% Add source term and divide by mass
dT = (dT + s(x))./m(x);
end