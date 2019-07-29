function config = ybranch_init()

config = struct;

config.C_OUTS = 2;
config.outs_names = {'lo';'hi'};

config.C_CONTS = 4;
config.cont_names = {'bot_out';'bot_in';'top_in';'top_out'};
config.signs = [1 -1 1 -1];

config.geom_script = char('ysplit_geom;');
config.adj_script = char('run_double_adjoint;');
config.fwd_script = char('run_double_forward;');

config.baseconts = yb_contours(7e-9);

config.tapers_fine = cell(config.C_CONTS,1);
for k = 1:config.C_CONTS
    config.tapers_fine{k} = [sin(linspace(0,pi/2,7)).^2 ones(1,size(config.baseconts{k},1)-14) cos(linspace(0,pi/2,7)).^2];
end

end
