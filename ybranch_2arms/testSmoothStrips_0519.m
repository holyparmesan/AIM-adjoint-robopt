
% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');
h = appopen('mode');
appevalscript(h,char('load("ybranch.lms");'));
appevalscript(h,char('addpath("../");'));

config = ybranch_init;
conts_fine = config.baseconts;

[delts_fine, norms_fine] = delts_norms(conts_fine{:});

N_ROUGHS = 5;
N_SMOOTHS = 5;
N_PER = 20;

ROUGH_AMPL = 10e-9;
ROUGH_CL = 40e-9;

SMOOTH_AMPL = 10e-9;
SMOOTH_CL = 400e-9;

scales = linspace(-1,1,N_PER);
fwd_flucts = cell(config.C_CONTS,1);
adj_flucts = cell(config.C_CONTS,1);
roughs_base = cell(config.C_CONTS,1);
smooths_base = cell(config.C_CONTS,1);

for k = 1:config.C_CONTS
    L_FINE = size(conts_fine{k},1);
    taper_rough = [sin(linspace(0,pi/2,7)).^2 ones(1,L_FINE - 14) cos(linspace(0,pi/2,7)).^2];
    roughs_base{k} = ROUGH_AMPL * randomRoughness(conts_fine{k},ROUGH_CL,N_ROUGHS) .* taper_rough;
    if or(k == 1, k == 4)
        taper_smooth = [sin(linspace(0,pi/2,15)).^2 ones(1,1000 - 30) cos(linspace(0,pi/2,15)).^2 zeros(1,L_FINE - 1000)];
        smooths_base{k} = SMOOTH_AMPL * randomRoughness(conts_fine{k},SMOOTH_CL,N_SMOOTHS) .* taper_smooth;
    else
        smooths_base{k} = zeros(N_SMOOTHS,L_FINE);
    end
    
    adj_flucts{k} = [zeros(1,L_FINE); roughs_base{k}];
    fwd_flucts{k} = zeros(N_ROUGHS*N_SMOOTHS*N_PER,L_FINE);
    for m = 1:N_ROUGHS
        for i = 1:N_PER
            fwd_flucts{k}(i+(m-1)*(N_PER*N_SMOOTHS):N_PER:i+m*(N_PER*N_SMOOTHS)-1,:) = roughs_base{k}(m,:) + smooths_base{k} * scales(i);
        end
    end
end

[hist_adj, hist_ssa] = run_adjes(h, config, adj_flucts);
hist_ssf = run_fwds(h, config, fwd_flucts);

save("ybranch_2arms/data/testSS_0519.mat");
