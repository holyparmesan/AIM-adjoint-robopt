function [y_bend, x_bend] = yb_out(t_l)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% bend
L = 10.5;
S = 2.4;
R = 4.5;
w = 0.5;
x_c = @(t) 3*(1-t).*(1-2*t).*t*R + t.^2.*(3-2*t)*L;
y_c = @(t) t.^2.*(3-2*t)*S;
dydx_c = @(t) 2*t.*(1-t)*S./ ((1-6.*t+6.*t.^2)*R+2*t.*(1-t)*L);
x_u = @(t) x_c(t) - w/2*sin(atan(dydx_c(t)));
y_u = @(t) y_c(t) + w/2*cos(atan(dydx_c(t)));

x_bend = x_u(t_l)+2;
y_bend = y_u(t_l)+0.35;
end