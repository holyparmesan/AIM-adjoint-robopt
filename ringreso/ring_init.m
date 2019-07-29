function config = ring_init(spacing)

config = struct;

config.C_OUTS = 1;
config.outs_names = {'hi'};

config.C_CONTS = 2;
config.cont_names = {'r_out';'r_in'};
config.signs = [1 -1];

config.geom_script = char('ring_geom;');
config.adj_script = char('run_badjoint;');
config.fwd_script = char('run_bforward;');

config.baseconts = basecont_ring(spacing);

end

