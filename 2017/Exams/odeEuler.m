function [ t,y ] = odeEuler( odefun, T, y0, dt )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n = ceil((T(2) - T(1))/dt);
t = linspace(T(1),T(2),n)';
y = zeros(length(t),length(y0));
y(1,:) = y0;
for i = 2:n
    y(i,:) = y(i-1,:) + dt*odefun(t(i-1),y(i-1,:))';
end
end

