
% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');
h = appopen('mode');
appevalscript(h,char('load("ybranch.lms");'));
appevalscript(h,char('addpath("../");'));

config = ybranch_init;
conts_fine = config.baseconts;
adjstr = run_adjoint(h,config,conts_fine);

N_CLS = 4;
N_LINES = 10;
N_PER = 20;

ROUGH_AMPL = 20e-9;
ROUGH_CLS = [20e-9, 50e-9, 100e-9, 400e-9];

N_MC = N_LINES * N_PER;
fwd_flucts = cell(config.C_CONTS,N_CLS);
scales = linspace(-1,1,N_PER);
flucts_base = cell(config.C_CONTS,N_CLS);

hists_ssf = cell(N_CLS,1);

for l = 1:N_CLS
    for k = 1:config.C_CONTS
        L_FINE = size(conts_fine{k},1);
        taper_fine = [sin(linspace(0,pi/2,7)).^2 ones(1,L_FINE - 14) cos(linspace(0,pi/2,7)).^2];
        flucts_base{k,l} = ROUGH_AMPL * randomRoughness(conts_fine{k},ROUGH_CLS(l),N_LINES) .* taper_fine;
        fwd_flucts{k,l} = zeros(N_MC,L_FINE);
        for i = 1:N_PER
            fwd_flucts{k,l}(i:N_PER:end,:) = flucts_base{k,l} * scales(i);
        end
    end

    save("testbk.mat");
    hists_ssf{l} = run_fwds(h, config, {fwd_flucts{:,l}});
end


save("ybranch_2arms/data/testSR_0520.mat");