function [surf_y, surf_x] = yb_in(t_l)
%UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% surf_x = surf_x(surf_x >= 2);
% surf_y = zeros(size(surf_x));
% 
% range_output = surf_x >= 12.5;
% surf_y(range_output) = 2.5;

% bend
L = 10.5;
S = 2.4;
R = 4.5;
w = 0.5;
x_c = @(t) 3*(1-t).*(1-2*t).*t*R + t.^2.*(3-2*t)*L;
y_c = @(t) t.^2.*(3-2*t)*S;
dydx_c = @(t) 2*t.*(1-t)*S./ ((1-6.*t+6.*t.^2)*R+2*t.*(1-t)*L);
x_l = @(t) x_c(t) + w/2*sin(atan(dydx_c(t)));
y_l = @(t) y_c(t) - w/2*cos(atan(dydx_c(t)));

% range_bend = surf_x < 12.5;
% bend_x = surf_x(range_bend) - 2;
% f = @(t) x_l(t) - bend_x;
% options = optimoptions('fsolve','Display','off');
% t_l = fsolve(f, bend_x/L, options);
surf_y = y_l(t_l) + 0.35;
surf_x = x_l(t_l) + 2;
end