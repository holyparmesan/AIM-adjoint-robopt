function [y_junc, x_junc] = yb_out(t_l)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% junction part
tape_x = linspace(0, 2, 13);
tape_y = [0.5 0.5 0.6 0.7 0.9 1.26 1.4 1.4 1.4 1.4 1.31 1.2 1.2]/2;
pp = spline(tape_x, tape_y);
x_junc = 2*t_l;
y_junc = ppval(pp, 2*t_l);

end