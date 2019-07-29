function adjstr = run_adjoint(h, config, rconts, varargin)
% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');
% rconts are the contours along which the adjoint is calculated.
% if varargin has content, it can be used to assist with cosining
% by providing a direction of travel which may not be normal to rconts.

% config is a struct which must have:
% C_OUTS and a cell outs_names{C_OUTS}
% C_CONTS and a cell cont_names{C_CONTS}
% geom_script and adj_script, which are char('')

for k = 1:config.C_CONTS
    appputvar(h,config.cont_names{k},rconts{k});
end

appevalscript(h,config.geom_script); % geo sensitive
appevalscript(h,config.adj_script); % geo sensitive

Edats = load("Epass.mat");
xs = appgetvar(h,char('xs'));
ys = appgetvar(h,char('ys'));
freqs = appgetvar(h,char('freqs'));
Edat_fwd = Edats.Edat_in; % appgetvar(h,char('Edat_in'));
idx_x = appgetvar(h,char('index_x'));
idx_y = appgetvar(h,char('index_y'));

N_FREQS = numel(freqs);

Ss = zeros(config.C_OUTS,N_FREQS);
adjSs = cell(config.C_CONTS,1);
for k = 1:config.C_CONTS
    adjSs{k} = zeros(config.C_OUTS,size(rconts{k},1),N_FREQS);
end

motionCosines = cell(config.C_CONTS,1);
if numel(varargin) > 0
    normsMove = varargin{1};
    [~, normsAdj] = delts_norms(rconts{:});
    for k = 1:config.C_CONTS
        motionCosines{k} = sum(normsMove{k}.' .* normsAdj{k}.',1);
    end
else
    for k = 1:config.C_CONTS
        motionCosines{k} = ones(1,size(rconts{k},1));
    end
end

for l = 1:config.C_OUTS
    Ss(l,:) = appgetvar(h,char(strcat('S_',config.outs_names{l}))).';
    Edat_a = getfield(Edats,strcat('Edat_',config.outs_names{l}));
%    Edat_a{l} = appgetvar(h,char(strcat('Edat_',config.outs_names{l})));
    Nmodes{l} = appgetvar(h,char(strcat('N_',config.outs_names{l})));
    adjust_a{l} = appgetvar(h,char(strcat('adjust_',config.outs_names{l}))).';

    adjraw = adjCalc(xs,ys,Edat_fwd,Edat_a,idx_x,idx_y,freqs,Nmodes{l},true,rconts{:});
    for k = 1:config.C_CONTS
        adjSs{k}(l,:,:) = config.signs(k) * motionCosines{k}.' .* adjraw{k} ./ adjust_a{l};
    end
end
clear Edat_fwd Edat_a;

adjstr.ss = Ss;
adjstr.adjss = adjSs;
adjstr.freqs = freqs;