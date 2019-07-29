function [output_log, mesh] = runSampsDC(profiles)

%% this script is for evaluating a directional coupler repeatedly

% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');

h = appopen('mode');
appevalscript(h,char('load("dircoupler.lms");'));

N_EVAL = size(profiles,1);

output_log = zeros(N_EVAL,10);

for i = 1:N_EVAL
    appputvar(h,char('bot_out'),squeeze(profiles(i,1,:)));
    appputvar(h,char('bot_in'),squeeze(profiles(i,2,:)));
    appputvar(h,char('top_in'),squeeze(profiles(i,3,:)));
    appputvar(h,char('top_out'),squeeze(profiles(i,4,:)));
    
    appevalscript(h,char('dircoupler_geom;'));
    appevalscript(h,char('run_varfdtd_s12;'));
    
    s12_data = appgetvar(h,char('S12'));
    mesh = appgetvar(h,char('m_a'));
    output_log(i,:) = s12_data;
end
    
appclose(h);

end