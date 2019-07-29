function [hist_adjss, hist_ss] = run_adjes(h, config, flucts)

conts_fine = config.baseconts;

[~, norms_fine] = delts_norms(conts_fine{:});

N_MC = size(flucts{1},1);
hist_adjss = cell(config.C_CONTS,1);
for i = 1:N_MC
    modconts = cell(config.C_CONTS,1);
    for k = 1:config.C_CONTS
        modconts{k} = conts_fine{k} + (flucts{k}(i,:).' .* norms_fine{k});
    end
    
    adjstr = run_adjoint(h,config,modconts,norms_fine);
    
    if i == 1
        N_FREQS = numel(adjstr.freqs);
        hist_ss = zeros(N_MC,config.C_OUTS,N_FREQS);
        
        for k = 1:config.C_CONTS
            L_FINE = size(flucts{k},2);
            hist_adjss{k} = zeros(N_MC,config.C_OUTS,L_FINE,N_FREQS);
        end
    end
    hist_ss(i,:,:) = adjstr.ss;
    for k = 1:config.C_CONTS
        hist_adjss{k}(i,:,:,:) = adjstr.adjss{k};
    end
    
    if mod(i,20) == 0
        save("databackup.mat");
    end
end

save("datafinal.mat");