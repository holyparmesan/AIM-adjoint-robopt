function adjFigMaker(conts_fine, adjss, signs, ssbase, freqind, varargin)
% varargin is for if you gotta index into adjss

N_DS = 10;

[delts, ~] = delts_norms(conts_fine{:});

cmp = diverging_map(linspace(0,1,100),[0.23, 0.299, 0.754],[0.706,0.016,0.150]);
max_limit = 4e12;
%cmp = jet;

 % conts_fine is cell(CONTS,1), internals are [L_FINE,2]
 % adjss is cell(CONTS,1), internals are [???, C_OUTS, L_FINE] where ???
 % may or may not exist
 % ssbase is [C_OUTS]. or so help me god.

f_main = figure;

subplot(2,1,1);
title("Magnitude (lower port)",'fontsize',12);
hold on;
view(0,90);
colormap(cmp);
caxis([-max_limit max_limit]);
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(findall(gca,'type','text'),'visible','on')
colorbar;

subplot(2,1,2);
title("Phase (lower port)",'fontsize',12);
hold on;
view(0,90);
colormap(cmp);
caxis([-max_limit max_limit]);
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(findall(gca,'type','text'),'visible','on')
colorbar;

f_diffs = figure;

subplot(2,1,1);
title("Power (total)",'fontsize',12);
hold on;
view(0,90);
colormap(cmp);
caxis([-max_limit max_limit]);
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(findall(gca,'type','text'),'visible','on')
colorbar;

subplot(2,1,2);
title("Power (imbalance)",'fontsize',12);
hold on;
view(0,90);
colormap(cmp);
caxis([-max_limit max_limit]);
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(findall(gca,'type','text'),'visible','on')
colorbar;
 
for k = 1:size(conts_fine,1)
    exes = conts_fine{k}(1:N_DS:end,1);
    whys = conts_fine{k}(1:N_DS:end,2);
    N_DSD = size(exes,1);
    
    if numel(varargin) > 0
        adj = signs(k) * shiftdim(adjss{k}(varargin{1},:,1:N_DS:end,freqind),1) ./ vecnorm(delts{k}(1:N_DS:end,:).');
    else
        adj = signs(k) * adjss{k}(:,1:N_DS:end,freqind) ./ vecnorm(delts{k}(1:N_DS:end,:).');
    end
    
    adj_re = real(adj(1,:).' ./ ssbase(1)) * abs(ssbase(1));
    adj_im = imag(adj(1,:).' ./ ssbase(1));
    adj_mag1 = 2 * real(adj(1,:) .* conj(ssbase(1))).';
    adj_mag2 = 2 * real(adj(2,:) .* conj(ssbase(2))).';
    
    figure(f_main);
    subplot(2,1,1);
    surf([exes,exes],[whys,whys],zeros(N_DSD,2),[adj_re,adj_re],'facecol','no','edgecol','interp','linew',4);
    subplot(2,1,2);
    surf([exes,exes],[whys,whys],zeros(N_DSD,2),[adj_im,adj_im],'facecol','no','edgecol','interp','linew',4);
    figure(f_diffs);
    subplot(2,1,1);
    surf([exes,exes],[whys,whys],zeros(N_DSD,2),[adj_mag1 + adj_mag2,adj_mag1 + adj_mag2],'facecol','no','edgecol','interp','linew',4);
    subplot(2,1,2);
    surf([exes,exes],[whys,whys],zeros(N_DSD,2),[adj_mag1 - adj_mag2,adj_mag1 - adj_mag2],'facecol','no','edgecol','interp','linew',4);
    
end

end