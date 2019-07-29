function hist_ss = run_fwds(h, config, flucts)

conts_fine = config.baseconts;
[~, norms_fine] = delts_norms(conts_fine{:});

N_MC = size(flucts{1},1);
modconts = cell(config.C_CONTS,1);

for i = 1:N_MC
    for k = 1:config.C_CONTS
        modconts{k} = conts_fine{k} + (flucts{k}(i,:).' .* norms_fine{k});
    end
    
    ss = run_forward(h,config,modconts);
    
    if i == 1
        N_FREQS = size(ss,2);
        hist_ss = zeros(N_MC,config.C_OUTS,N_FREQS);
    end
    
    hist_ss(i,:,:) = ss;
    
    if mod(i,20) == 0
        save("databackup.mat");
    end
end

save("datafinal.mat");