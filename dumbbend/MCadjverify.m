function ret = MCadjverify(testcontours, N_MC, varargin)
% testcontours: N_DESNS by L_DESIGN by 2,
% If you want to only test 1 contour, use shiftdim(cont,-1) before calling
% this.

% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');

ret = struct;

h = appopen('mode');
appevalscript(h,char('load("ybranch.lms");'));
appevalscript(h,char('addpath("../");'));

config = ybranch_init();

appputvar(h,char('bot_in'),config.baseconts{2});
appputvar(h,char('top_in'),config.baseconts{3});

N_DESNS = size(testcontours,1);
L_DESIGN = size(testcontours,2);
taper_desn = [sin(linspace(0,pi/2,7)).^2 ones(1,L_DESIGN-14) cos(linspace(0,pi/2,7)).^2];

L_INTERP = 5;
L_FINE = (L_DESIGN-1) * L_INTERP + 1;
taper_fine = interp1(linspace(0,1,L_DESIGN),taper_desn,linspace(0,1,L_FINE),'spline');

conts_fine_lo = zeros(N_DESNS, L_FINE, 2);
for i = 1:N_DESNS
    conts_fine_lo(i,:,1) = interp1(linspace(0,1,L_DESIGN),testcontours(i,:,1),linspace(0,1,L_FINE),'spline');
    conts_fine_lo(i,:,2) = interp1(linspace(0,1,L_DESIGN),testcontours(i,:,2),linspace(0,1,L_FINE),'spline');
end

delts_fine = zeros(size(conts_fine_lo));
delts_fine(:,1,:) = conts_fine_lo(:,2,:) - conts_fine_lo(:,1,:);
delts_fine(:,2:end-1,:) = conts_fine_lo(:,3:end,:) - conts_fine_lo(:,1:end-2,:);
delts_fine(:,end,:) = conts_fine_lo(:,end,:) - conts_fine_lo(:,end-1,:);

norms_mags = sqrt(delts_fine(:,:,1).^2 + delts_fine(:,:,2).^2);

norms_fine_lo = zeros(size(conts_fine_lo));
norms_fine_lo(:,:,1) = delts_fine(:,:,2) ./ norms_mags;
norms_fine_lo(:,:,2) = -1*delts_fine(:,:,1) ./ norms_mags;

conts_fine_hi = zeros(size(conts_fine_lo));
conts_fine_hi(:,:,1) = conts_fine_lo(:,:,1);
conts_fine_hi(:,:,2) = -1*conts_fine_lo(:,:,2);

norms_fine_hi = zeros(size(conts_fine_lo));
norms_fine_hi(:,:,1) = norms_fine_lo(:,:,1);
norms_fine_hi(:,:,2) = -1*norms_fine_lo(:,:,2);

ROUGH_AMPL = 10e-9;
ROUGH_CL = 50e-9;

if nargin>2
    ret.scales = varargin{1};
    N_REPEAT = numel(ret.scales);
    ret.flucts_base = ROUGH_AMPL * randomRoughness(squeeze(conts_fine_lo(1,:,:)),ROUGH_CL,N_MC*2); % probabilities according to base - shouldn't matter too much
    N_MC = N_MC * N_REPEAT;
    flucts = zeros((N_MC+1)*2,L_FINE);
    for i = 1:N_REPEAT
        flucts(1+2*i:2*N_REPEAT:end,:) = ret.flucts_base(1:2:end,:) * ret.scales(i);
        flucts(2+2*i:2*N_REPEAT:end,:) = ret.flucts_base(2:2:end,:) * ret.scales(i);
    end
else
    flucts = ROUGH_AMPL * randomRoughness(squeeze(conts_fine_lo(1,:,:)),ROUGH_CL,(N_MC+1)*2); % probabilities according to base - shouldn't matter too much
    flucts(1:2,:) = zeros(2,L_FINE);
    N_REPEAT = 1;
end
    
ret.sslo = zeros(N_MC+1,N_DESNS);
ret.sshi = zeros(N_MC+1,N_DESNS);

for i = 1:N_MC+1
    for j = 1:N_DESNS        
        modlo = conts_fine_lo(j,:,:) + (taper_fine .* flucts(2*i-1,:) .* norms_fine_lo(j,:,:));
        modhi = conts_fine_hi(j,:,:) + (taper_fine .* flucts(2*i,:) .* norms_fine_hi(j,:,:));
        
        appputvar(h,config.cont_names{1},squeeze(modlo));
        appputvar(h,config.cont_names{2},squeeze(modhi));
        
        appevalscript(h,char('ysplit_geom;'));
        
        appevalscript(h,char('run_double_forward;'));
        ret.sslo(i,j,:) = appgetvar(h,char('S_lo'));
        ret.sshi(i,j,:) = appgetvar(h,char('S_hi'));          
    end
end

ret.imbals = (abs(ret.sslo).^2 - abs(ret.sshi).^2);
ret.losses = 1 - (abs(ret.sslo).^2 + abs(ret.sshi).^2);
ret.flucts = flucts;
ret.conts_fine_lo = conts_fine_lo;

save("databackupMC.mat");

end