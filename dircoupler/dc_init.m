function config = dc_init(L_FINE, inside, outside)
% inside = 0.1e-6, outside = 0.6e-6 is standard

config = struct;

config.C_OUTS = 2;
config.outs_names = {'lo';'hi'};

config.C_CONTS = 4;
config.cont_names = {'bot_out';'bot_in';'top_in';'top_out'};
config.signs = [1 -1 1 -1];

config.N_FREQS = 9; % careful...

config.geom_script = char('dircoupler_geom_new;');
config.adj_script = char('run_double_adjoint;');
config.fwd_script = char('run_double_forward;');

x0s = [linspace(-1.5e-6,1.5e-6,L_FINE);zeros(1,L_FINE)].';
y1s = [zeros(1,L_FINE);ones(1,L_FINE)].';
config.baseconts = {x0s - y1s*outside;x0s - y1s*inside;x0s + y1s*inside;x0s + y1s*outside};

end