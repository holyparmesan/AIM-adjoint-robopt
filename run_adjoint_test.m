function adjstr = run_adjoint_test(h, config, rconts, travelNorms)
% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');
% rconts are the contours along which the adjoint is calculated.
% if varargin has content, it can be used to specify alternate directions
% for adjoint calculation.

% config is a struct which must have:
% C_OUTS and a cell outs_names{C_OUTS}
% C_CONTS and a cell cont_names{C_CONTS}
% geom_script and adj_script, which are char('')

% The point of this routine is to test adjCalc vs. adjCalcFIT, and also
% cosining vs. no cosining.

for k = 1:config.C_CONTS
    appputvar(h,config.cont_names{k},rconts{k});
end

appevalscript(h,config.geom_script); % geo sensitive
appevalscript(h,config.adj_script); % geo sensitive

xs = appgetvar(h,char('xs'));
ys = appgetvar(h,char('ys'));
freqs = appgetvar(h,char('freqs'));
Edat_fwd = appgetvar(h,char('Edat_in'));
idx_x = appgetvar(h,char('index_x'));
idx_y = appgetvar(h,char('index_y'));

N_FREQS = numel(freqs);

Ss = zeros(config.C_OUTS,N_FREQS);
adjSs_int = cell(config.C_CONTS);
adjSs_fit = cell(config.C_CONTS);
for k = 1:config.C_CONTS
    adjSs_int{k} = zeros(config.C_OUTS,size(rconts{k},1),N_FREQS);
    adjSs_fit{k} = zeros(config.C_OUTS,size(rconts{k},1),N_FREQS);
end

for l = 1:config.C_OUTS
    Ss(l,:) = appgetvar(h,char(strcat('S_',config.outs_names{l}))).';
    Edat_a{l} = appgetvar(h,char(strcat('Edat_',config.outs_names{l})));
    Nmodes{l} = appgetvar(h,char(strcat('N_',config.outs_names{l})));
    adjust_a{l} = appgetvar(h,char(strcat('adjust_',config.outs_names{l}))).';

    adjraw_int = adjCalc(xs,ys,Edat_fwd,Edat_a{l},idx_x,idx_y,freqs,Nmodes{l},true,rconts{:});
    adjraw_fit = adjCalcFIT(xs,ys,Edat_fwd,Edat_a{l},idx_x,idx_y,freqs,Nmodes{l},true,rconts{:});
    for k = 1:config.C_CONTS
        adjSs_int{k}(l,:,:) = adjraw_int{k} ./ adjust_a{l};
        adjSs_fit{k}(l,:,:) = adjraw_fit{k} ./ adjust_a{l};
        % this is the place where you would multiply by cosines{k}
    end
           
    clear Edat_in Edat_a;
end

adjstr.cosines = cell(config.C_CONTS,1);
[~, normsAdj] = delts_norms(rconts{:});
for k = 1:config.C_CONTS
    adjstr.cosines{k} = sum(travelNorms{k}.' .* normsAdj{k}.',1);
end

adjstr.ss = Ss;
adjstr.adjss_int = adjSs_int;
adjstr.adjss_fit = adjSs_fit;
adjstr.freqs = freqs;

