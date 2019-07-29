% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');
h = appopen('mode');
appevalscript(h,char('load("dircoupler_narrow.lms");'));
appevalscript(h,char('addpath("../");'));

config = dc_init(200, 0.075e-6, 0.575e-6);
conts_fine = config.baseconts;

adjstr = run_adjoint(h,config,conts_fine);
ssbase = adjstr.ss;
adjss = adjstr.adjss;

[delts_fine, norms_fine] = delts_norms(conts_fine{:});

N_CLS = 3;
N_LINES = 10;
N_PER = 20;

ROUGH_AMPL = 20e-9;
ROUGH_CLS = [20e-9, 50e-9, 100e-9];

N_MC = N_LINES * N_PER;
flucts = cell(config.C_CONTS,N_CLS);
scales = linspace(-1,1,N_PER);
flucts_base = cell(config.C_CONTS,N_CLS);
hists_ss = cell(N_CLS,1);

for l = 1:N_CLS
    save("testbk.mat");
    for k = 1:config.C_CONTS
        L_FINE = size(conts_fine{k},1);
        taper_fine = taper(L_FINE, 7);
        flucts_base{k,l} = ROUGH_AMPL * randomRoughness(conts_fine{k},ROUGH_CLS(l),N_LINES) .* taper_fine;
        flucts{k,l} = zeros(N_MC,L_FINE);
        for i = 1:N_PER
            flucts{k,l}(i:N_PER:end,:) = flucts_base{k,l} * scales(i);
        end
    end

    hists_ss{l} = run_fwds(h, config, {flucts{:,l}});
end

save("test0509.mat");

function tap = taper(L, l_taper)
    tap = [sin(linspace(0,pi/2,l_taper)).^2 ones(1,L - 2*l_taper) cos(linspace(0,pi/2,l_taper)).^2];
end