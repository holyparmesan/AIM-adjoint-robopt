function baseconts = basecont_gen(t_l)
% the understanding is t_l goes from 0 to 1.

% bend geometry
L = 5.0;
H = 3.0;
w = 0.5;

x_c = @(t) L * t;
y_c = @(t) H * sin(pi*t/2).^2;
dydx_c = @(t) sin(pi*t) * (pi * H) / (2 * L);
dx_l = @(t) w/2*sin(atan(dydx_c(t)));
dy_l = @(t) -w/2*cos(atan(dydx_c(t)));

y_lo = y_c(t_l) + dy_l(t_l) - H/2;
y_hi = y_c(t_l) - dy_l(t_l) - H/2;
x_lo = x_c(t_l) + dx_l(t_l) - L/2;
x_hi = x_c(t_l) - dx_l(t_l) - L/2;

baseconts = cell(2,1);
baseconts{1} = 1e-6 * [x_lo ; y_lo].';
baseconts{2} = 1e-6 * [x_hi ; y_hi].';

end