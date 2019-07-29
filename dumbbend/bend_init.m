function config = bend_init(t_l)

config = struct;

config.C_OUTS = 1;
config.outs_names = {'out'};

config.C_CONTS = 2;
config.cont_names = {'bot';'top'};
config.signs = [1 -1];

config.geom_script = char('bend_geom;');
config.adj_script = char('run_adjoint;');
config.fwd_script = char('run_forward;');

config.baseconts = basecont_gen(t_l);

end

