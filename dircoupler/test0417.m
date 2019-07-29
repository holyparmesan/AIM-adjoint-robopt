
% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');
h = appopen('mode');
appevalscript(h,char('load("dircoupler_wide.lms");'));
appevalscript(h,char('addpath("../");'));

L_FINE = 100;
config = dc_init(L_FINE, 0.15e-6, 0.65e-6);
conts_fine = config.baseconts;

[delts_fine, norms_fine] = delts_norms(conts_fine{:});

N_LINES = 20;
N_PER = 20;

widths = 2 * floor(linspace(6,L_FINE-14,N_LINES).'/2);
tapers_fine = zeros(N_LINES,L_FINE);
for i = 1:N_LINES
    n = widths(i);
    tapers_fine(i,:) = [zeros(1, (L_FINE - (14+n))/2) sin(linspace(0,pi/2,7)).^2 ones(1,n) cos(linspace(0,pi/2,7)).^2 zeros(1, (L_FINE - (14+n))/2)];
end

ROUGH_AMPL = 10e-9;
ROUGH_CL = 50e-9;

N_MC = N_LINES * N_PER;
flucts = cell(config.C_CONTS,1);
scales = linspace(-1,1,N_PER);
flucts_base = cell(config.C_CONTS,1);

for k = 1:config.C_CONTS
    flucts_base{k} = ROUGH_AMPL * randomRoughness(conts_fine{k},ROUGH_CL,N_LINES) .* tapers_fine;
    flucts{k} = zeros(1+N_MC,L_FINE);
    for i = 1:N_PER
        flucts{k}(1+i:N_PER:end,:) = flucts_base{k} * scales(i);
    end
end

N_FREQS = 0;
hist_adjss = cell(config.C_CONTS,1);

for i = 1:N_MC+1
    modconts = cell(config.C_CONTS,1);
    for k = 1:config.C_CONTS
        modconts{k} = conts_fine{k} + (flucts{k}(i,:).' .* norms_fine{k});
    end
    
    adjstr = run_adjoint(h,config,modconts);  
    if i == 1
        N_FREQS = numel(adjstr.freqs);
        hist_ss = zeros(N_MC+1,config.C_OUTS,N_FREQS);
        for k = 1:config.C_CONTS
            hist_adjss{k} = zeros(N_MC+1,config.C_OUTS,L_FINE,N_FREQS);
        end
    end
    hist_ss(i,:,:) = adjstr.ss;
    for k = 1:config.C_CONTS
        hist_adjss{k}(i,:,:,:) = adjstr.adjss{k};
    end
    
    if mod(i,20) == 0
        save("databackup.mat");
    end
end

save("datafinal.mat");