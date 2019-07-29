
% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');
h = appopen('mode');
appevalscript(h,char('load("dumbtest.lms");'));
appevalscript(h,char('addpath("../");'));

config = bend_init(linspace(0,1,75));

L_DESIGN = 75;
taper_desn = [sin(linspace(0,pi/2,7)).^2 ones(1,L_DESIGN-14) cos(linspace(0,pi/2,7)).^2];

L_INTERP = 5;
L_FINE = (L_DESIGN-1) * L_INTERP + 1;
taper_fine = interp1(linspace(0,1,L_DESIGN),taper_desn,linspace(0,1,L_FINE),'spline');

N_LINES = 10;
N_PER = 30;

conts_fine = cell(2,1);
conts_fine{1} = zeros(L_FINE, 2);
conts_fine{2} = zeros(L_FINE, 2);
for i = 1:2
    conts_fine{i}(:,1) = interp1(linspace(0,1,L_DESIGN),config.baseconts{i}(:,1),linspace(0,1,L_FINE),'spline');
    conts_fine{i}(:,2) = interp1(linspace(0,1,L_DESIGN),config.baseconts{i}(:,2),linspace(0,1,L_FINE),'spline');
end

[delts_fine, norms_fine] = delts_norms(conts_fine{:});

ROUGH_AMPL = 10e-9;
ROUGH_CL = 50e-9;

N_MC = N_LINES * N_PER;
flucts = zeros(1+N_MC,2,L_FINE);
flucts_base_lo = ROUGH_AMPL * randomRoughness(conts_fine{1},ROUGH_CL,N_LINES);
flucts_base_hi = ROUGH_AMPL * randomRoughness(conts_fine{2},ROUGH_CL,N_LINES);
scales = linspace(-1,1,N_PER);

for i = 1:N_PER
    flucts(1+i:N_PER:end,1,:) = flucts_base_lo * scales(i);
    flucts(1+i:N_PER:end,2,:) = flucts_base_hi * scales(i);
end

N_FREQS = 0;
hist_adjss = cell(config.C_CONTS,1);

for i = 1:N_MC+1
    modconts = cell(config.C_CONTS,1);
    for k = 1:config.C_CONTS
        modconts{k} = conts_fine{k} + (taper_fine.' .* squeeze(flucts(i,k,:)) .* norms_fine{k});
    end
    
    adjstr = run_adjoint(h,config,modconts);  
    if i == 1
        N_FREQS = numel(adjstr.freqs);
        hist_ss = zeros(N_MC+1,config.C_OUTS,N_FREQS);
        for k = 1:config.C_CONTS
            hist_adjss{k} = zeros(N_MC+1,config.C_OUTS,L_FINE,N_FREQS);
        end
    end
    hist_ss(i,:) = adjstr.ss;
    for k = 1:config.C_CONTS
        hist_adjss{k}(i,:,:,:) = adjstr.adjss{k};
    end
    
    if mod(i,20) == 0
        save("databackup.mat");
    end
end

save("databackup.mat");