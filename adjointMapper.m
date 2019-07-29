function [adj_maps, lambs_s] = adjointMapper(h, varargin)
% each varargin is a contour: an N x 2 list of vertices

appevalscript(h,char('run_adjoint;'));

lambs_s = appgetvar(h,char('lambs_s'));
lambs_mon = appgetvar(h,char('lambs_mon'));

if (size(lambs_s,1) ~= size(lambs_mon,1)) || norm(lambs_s - lambs_mon) > 1e-10
    disp("Warning - mismatched port frequencies and monitor frequencies!");
end

xs = appgetvar(h,char('xs'));
ys = appgetvar(h,char('ys'));

E12 = appgetvar(h,char('Edat12'));
E21 = appgetvar(h,char('Edat21'));

ix = appgetvar(h,char('index_x'));
iy = appgetvar(h,char('index_y'));
iz = appgetvar(h,char('index_z'));

n_lo = 1.444;
n_hi = 3.4766;

adj_maps = adjCalc(xs, ys, E12, E21, ix, iy, iz, n_lo, n_hi, varargin{:});

end