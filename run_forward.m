function Ss = run_forward(h, config, rconts)
% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');
% rconts are the contours along which the shape is altered.

% config is a struct which must have:
% C_OUTS and a cell outs_names{C_OUTS}
% C_CONTS and a cell cont_names{C_CONTS}
% geom_script and fwd_script, which are char('')

for k = 1:config.C_CONTS
    appputvar(h,config.cont_names{k},rconts{k});
end

appevalscript(h,config.geom_script); % geo sensitive
appevalscript(h,config.fwd_script); % geo sensitive

freqs = appgetvar(h,char('freqs'));

Ss = zeros(config.C_OUTS,numel(freqs));

for l = 1:config.C_OUTS
    Ss(l,:) = appgetvar(h,char(strcat('S_',config.outs_names{l}))).';
end