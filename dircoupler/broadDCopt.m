% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');

N_ITERS = 50;
N_CONTOUR = 150;
N_FREQ = 9;

alpha = 8e-2 * 3e-8 * 1e-6 * 1e-6; % proportionality constant from adjoint field (1/m^2) to step size (m)
taper = [sin(linspace(0,pi/2,15)).^2 ones(1,N_CONTOUR-30) cos(linspace(0,pi/2,15)).^2];

contours = cell(4,1);
gradients = cell(4,1);
baseheights = [-0.6e-6, -0.1e-6, 0.1e-6, 0.6e-6];
signs = [1, -1, 1, -1];

hist_poses = cell(4,1);
hist_mapsl = cell(4,1);
hist_mapsh = cell(4,1);
hist_grads = cell(4,1);
hist_Shi = zeros(N_ITERS,N_FREQ);
hist_Slo = zeros(N_ITERS,N_FREQ);

for k = 1:4
   hist_poses{k} = zeros(N_ITERS,N_CONTOUR,2);
   hist_mapsl{k} = zeros(N_ITERS,N_CONTOUR,N_FREQ);
   hist_mapsh{k} = zeros(N_ITERS,N_CONTOUR,N_FREQ);
   hist_grads{k} = zeros(N_ITERS,N_CONTOUR);
   contours{k} = [linspace(-1.5e-6,1.5e-6,N_CONTOUR).' baseheights(k)*ones(N_CONTOUR,1)];
end

h = appopen('mode');
appevalscript(h,char('load("dircoupler.lms");'));
appevalscript(h,char('addpath("../");'));

for i = 1:N_ITERS
    for k = 1:4
        appputvar(h,char(strcat('contours',num2str(k))),contours{k});
    end
    appevalscript(h,char('dircoupler_geom_new;'));
    appevalscript(h,char('run_double_adjoint;'));

    xs = appgetvar(h,char('xs'));
    ys = appgetvar(h,char('ys'));
    freqs = appgetvar(h,char('freqs'));

    S_lo = appgetvar(h,char('S_lo')).';
    S_hi = appgetvar(h,char('S_hi')).';
    hist_Slo(i,:) = S_lo;
    hist_Shi(i,:) = S_hi;

    Edat_fwd = appgetvar(h,char('Edat_in'));
    Edat_ahi = appgetvar(h,char('Edat_hi'));
    Edat_alo = appgetvar(h,char('Edat_lo'));
    idx_x = appgetvar(h,char('index_x'));
    idx_y = appgetvar(h,char('index_y'));

    Nmode_lo = appgetvar(h,char('N_lo'));
    Nmode_hi = appgetvar(h,char('N_hi'));
    adjust_alo = appgetvar(h,char('adjust_lo')).';
    adjust_ahi = appgetvar(h,char('adjust_hi')).';

    adjmaps_lo_raw = adjCalc(xs,ys,Edat_fwd,Edat_alo,idx_x,idx_y,freqs,Nmode_lo,false,contours{:});
    adjmaps_hi_raw = adjCalc(xs,ys,Edat_fwd,Edat_ahi,idx_x,idx_y,freqs,Nmode_hi,false,contours{:});

    Tlo = abs(S_lo).^2;
    Thi = abs(S_hi).^2;
	% error function = (Tlo + Thi - 1)^2 + imbal_weight * (Tlo - Thi)^2
	% so gradient = 2(Tlo + Thi - 1)dT +/- imbal_weight * 2(Tlo - Thi)dT
    imbal_weight = 0.2;
    dTlo_coeff = 2 * (Tlo + Thi - 1 + imbal_weight * (Tlo - Thi));
    dThi_coeff = 2 * (Tlo + Thi - 1 - imbal_weight * (Tlo - Thi));

    for k = 1:4
        hist_poses{k}(i,:,:) = contours{k};
        hist_mapsl{k}(i,:,:) = adjmaps_lo_raw{k} ./ adjust_alo;
        hist_mapsh{k}(i,:,:) = adjmaps_hi_raw{k} ./ adjust_ahi;

        % T = ss*, so dT = 2 re{s' s*}
        gradients{k} = 2 * real(adjmaps_lo_raw{k} ./ adjust_alo .* conj(S_lo)) .* dTlo_coeff ...
            + 2 * real(adjmaps_hi_raw{k} ./ adjust_ahi .* conj(S_hi)) .* dThi_coeff;
        hist_grads{k}(i,:) = sum(gradients{k},2);
    end
    
    for k = 1:4
        % this is a hack, it should be perpendicular to the delta, not vertical. if we get serious curves this is going to break down.
        % contours{k}(:,2) = contours{k}(:,2) + signs(k) * alpha * taper.' .* sum(gradients{k},2);
        
        % this is a symmetric hack
        contours{k}(:,2) = contours{k}(:,2) + signs(k) * alpha * taper.' .* (sum(gradients{k},2) + sum(gradients{5-k},2)) * 0.5;
    end
    
end
